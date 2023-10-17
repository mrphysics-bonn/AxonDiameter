""" This function adds a reference scan to a sequence for calculation of sensitivity maps and B0 maps

    seq: Existing PyPulseq sequence
    prot: ISMRMRD metadata file
    system: MR system parameters
    params: Refscan parameters
    grads_off: Turn gradients off (only for ECC simulation)
"""
import math
import numpy as np

import ismrmrd

from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts

import pulseq_helper as ph

def gre_refscan_B0(seq, prot=None, system=Opts(), params=None, grads_off=False):

    # limit slew rate and maximum amplitude to avoid stimulation
    save_slew = system.max_slew
    save_grad = system.max_grad
    system.max_slew = 140 * system.gamma
    system.max_grad = 50e-3 * system.gamma

    # default parameters
    params_def = {"TE": [2.04e-3, 4.08e-3], "fov":210e-3, "res":2e-3, "ref_lines": 0, "flip_angle":12, "rf_dur":0.8e-3, 
                    "tbp": 2, "slices":1, "slice_res":2e-3, "dist_fac":0, "readout_bw": 1200, 
                    'separate_tr': False, 'center_out': False}
    for key in params:
        params_def[key] = params[key] # insert incoming parameters
    params = params_def
    
    # RF
    rf, gz, gz_reph, rf_del = make_sinc_pulse(flip_angle=params["flip_angle"] * math.pi / 180, duration=params["rf_dur"], slice_thickness=params["slice_res"],
                                apodization=0.5, time_bw_product=params["tbp"], system=system, return_gz=True, return_delay=True)

    # Calculate readout gradient and ADC parameters
    delta_k = 1 / params["fov"]
    if params["ref_lines"] > 0:
        Ny = int(params["ref_lines"])
    else:
        Ny = int(params["fov"]/params["res"]+0.5)
    if Ny % 2 != 0: 
        Ny -= 1 # even matrix size (+matrix size should not be bigger than imaging matrix)
    Nx = Ny
    samples = 2*Nx # 2x oversampling
    gx_flat_time_us = int(1e6/params["readout_bw"]) # readout_bw is in Hz/Px
    dwelltime = ph.trunc_to_raster(1e-6*gx_flat_time_us / samples, decimals=7)
    gx_flat_time = round(dwelltime*samples, 5)
    if (1e5*gx_flat_time %2 == 1):
        gx_flat_time += 10e-6 # even flat time
    diff_flat_adc = gx_flat_time - (dwelltime*samples)

    # Gradients
    gx_flat_area = Nx * delta_k * (gx_flat_time / (dwelltime*samples)) # compensate for longer flat time than ADC
    gx = make_trapezoid(channel='x', flat_area=gx_flat_area, flat_time=gx_flat_time, system=system)
    gx_mid = make_trapezoid(channel='x', area=-gx.area, system=system)
    gx_pre = make_trapezoid(channel='x', area=-gx.area / 2, system=system)
    gx_pre.delay = calc_duration(gz_reph) - calc_duration(gx_pre)
    phase_areas = (np.arange(Ny) - Ny / 2) * delta_k

    # spoilers
    gx_spoil = make_trapezoid(channel='x', area=2.5 * Nx * delta_k, system=system)
    gx_end = ph.merge_ramps([gx, gx_spoil], system=system) # merge x-spoiler with read gradient
    delay_end = calc_duration(gx)
    gz_spoil = make_trapezoid(channel='z', area=2/params['slice_res'], system=system)

    # calculate TE delays
    max_gy_pre = make_trapezoid(channel='y', area=max(abs(phase_areas)), system=system)
    gy_pre_dur = calc_duration(max_gy_pre)
    min_TE = np.ceil((gz.fall_time + gz.flat_time / 2 + calc_duration(gx_pre, max_gy_pre, gz_reph) + calc_duration(gx) / 2) / seq.grad_raster_time) * seq.grad_raster_time
    
    params["TE"] = sorted(params["TE"])
    n_TE = len(params["TE"])
    d_TE = np.asarray(params["TE"][1:]) - np.asarray(params["TE"][:-1])
    delay_TE = [params["TE"][0] - min_TE]
    for n_te, d_te in enumerate(d_TE):
        if not params['separate_tr']:
            delay_TE.append(d_te - calc_duration(gx) - calc_duration(gx_mid))
        else:
            delay_TE.append(params["TE"][n_te+1] - min_TE)
    if min(delay_TE) < 0:
        raise ValueError(f"TE is too small by {1e3*abs(min(delay_TE))} ms. Increase readout bandwidth.")

    # ADC 
    adc = make_adc(num_samples=samples, dwell=dwelltime, delay=gx.rise_time+diff_flat_adc/2, system=system)

    # RF spoiling
    rf_spoiling_inc = 117
    rf_phase = 0
    rf_inc = 0

    # build sequence
    if params['separate_tr']:
        TR = params['slices'] * (rf_del.delay + max(calc_duration(gx_pre,gz_reph),gy_pre_dur) + delay_TE[-1] 
             + calc_duration(gx) + max(gy_pre_dur,calc_duration(gz_spoil)))
    else:
        TR = params['slices'] * (rf_del.delay + max(calc_duration(gx_pre,gz_reph),gy_pre_dur) + delay_TE[0] 
             + (n_TE-1)*(calc_duration(gx) + calc_duration(gx_mid)) + np.sum(delay_TE[1:]) + max(gy_pre_dur+delay_end,calc_duration(gx_end)))

    print(f"Refscan volume TR: {TR*1e3:.3f} ms")

   # turn off gradients, if selected
    if grads_off:
        gz.amplitude = 0
        gz_reph.amplitude = 0
        gx_pre.amplitude = 0
        gx_spoil.amplitude = 0
        gz_spoil.amplitude = 0
        gx.amplitude = 0
        gx_mid.amplitude = 0
        gx_end.waveform = np.zeros_like(gx_end.waveform)
        phase_areas *= 0
    
    # center out readout
    if params['center_out']:
        phs_ix = Ny // 2
        prepscans = min(int(3/TR),25) # more dummy scans as we start in kspace center
    else:
        phs_ix = 0
        prepscans = min(int(1/TR),25)

    print(f"Refscan prepscans: {prepscans} ")

    if not params['separate_tr']:

        # imaging scans
        for i in range(-prepscans, Ny):
            ix = i if i >= 0 else None
            if ix is not None:
                if params['center_out']:
                    phs_ix += (-1)**(ix%2) * ix
                else:
                    phs_ix = ix
            if params["slices"]%2 == 1:
                slc = 0
            else:
                slc = 1

            # RF spoiling
            rf.phase_offset = rf_phase / 180 * np.pi
            adc.phase_offset = rf_phase / 180 * np.pi
            rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
            rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]

            for s in range(params["slices"]):
                if s==int(params["slices"]/2+0.5): 
                    if params["slices"]%2 == 1:
                        slc = 1
                    else:
                        slc = 0
                rf.freq_offset = gz.amplitude * params["slice_res"] * (slc - (params["slices"] - 1) / 2) * (1+params["dist_fac"]*1e-2)

                seq.add_block(rf, gz, rf_del)
                gy_pre = make_trapezoid(channel='y', area=phase_areas[phs_ix], duration=gy_pre_dur, system=system)
                seq.add_block(gx_pre, gy_pre, gz_reph)
                seq.add_block(make_delay(delay_TE[0]))
                for k in range(n_TE-1):
                    if i < 0:
                        seq.add_block(gx)
                    else:
                        seq.add_block(gx, adc)
                    seq.add_block(gx_mid)
                    seq.add_block(make_delay(delay_TE[k+1]))
                gy_pre.amplitude = -gy_pre.amplitude
                gy_pre.delay = delay_end
                if i < 0:
                    seq.add_block(gx_end, gy_pre)
                else:
                    seq.add_block(gx_end, adc, gy_pre)

                if prot is not None and i >= 0:
                    for k in range(n_TE):
                        acq = ismrmrd.Acquisition()
                        acq.idx.kspace_encode_step_1 = phs_ix
                        acq.idx.kspace_encode_step_2 = 0 # only 2D atm
                        acq.idx.slice = slc
                        acq.idx.contrast = k
                        acq.setFlag(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION)
                        if i == Ny-1 and k == n_TE-1:
                            acq.setFlag(ismrmrd.ACQ_LAST_IN_SLICE)
                        prot.append_acquisition(acq)
                    
                slc += 2 # interleaved
    else:

        # imaging scans
        for i in range(-prepscans, Ny):
            ix = i if i >= 0 else None
            if ix is not None:
                if params['center_out']:
                    phs_ix += (-1)**(ix%2) * ix
                else:
                    phs_ix = ix

            for k in range(n_TE):
                if params["slices"]%2 == 1:
                    slc = 0
                else:
                    slc = 1

                # RF spoiling
                rf.phase_offset = rf_phase / 180 * np.pi
                adc.phase_offset = rf_phase / 180 * np.pi
                rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
                rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]

                for s in range(params["slices"]):
                    if s==int(params["slices"]/2+0.5): 
                        if params["slices"]%2 == 1:
                            slc = 1
                        else:
                            slc = 0
                    rf.freq_offset = gz.amplitude * params["slice_res"] * (slc - (params["slices"] - 1) / 2) * (1+params["dist_fac"]*1e-2)

                    seq.add_block(rf, gz, rf_del)
                    gy_pre = make_trapezoid(channel='y', area=phase_areas[phs_ix], duration=gy_pre_dur, system=system)
                    seq.add_block(gx_pre, gy_pre, gz_reph)
                    seq.add_block(make_delay(delay_TE[k]))
                    if i < 0:
                        seq.add_block(gx)
                    else:
                        seq.add_block(gx, adc)
                    if k<n_TE-1:
                        seq.add_block(make_delay(d=params["TE"][-1]-params["TE"][k])) # account for longer 2nd echo, spoiler always in same position
                    gy_pre.amplitude = -gy_pre.amplitude
                    seq.add_block(gy_pre, gz_spoil)

                    if prot is not None and i>=0:
                        acq = ismrmrd.Acquisition()
                        acq.idx.kspace_encode_step_1 = phs_ix
                        acq.idx.kspace_encode_step_2 = 0 # only 2D atm
                        acq.idx.slice = slc
                        acq.idx.contrast = k
                        acq.setFlag(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION)
                        if i == Ny-1 and k == n_TE-1:
                            acq.setFlag(ismrmrd.ACQ_LAST_IN_SLICE)
                        prot.append_acquisition(acq)
                    
                    slc += 2
    
    delay_end = make_delay(d=3) # 3s delay after reference scan to allow for relaxation
    seq.add_block(delay_end)
    system.max_slew = save_slew
    system.max_grad = save_grad

    # use array to save echo times in protocol
    if prot is not None and n_TE>1:
        prot.append_array("echo_times", np.asarray(params["TE"]))
