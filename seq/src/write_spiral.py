# Spiral Pulseq Sequence

#%%

import numpy as np
import ismrmrd
import os
import datetime
from copy import copy

from pypulseq.make_arbitrary_grad import make_arbitrary_grad
from pypulseq.Sequence.sequence import Sequence
from pypulseq.make_adc import make_adc
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_gauss_pulse import make_gauss_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.make_delay import make_delay
from pypulseq.make_digital_output_pulse import make_digital_output_pulse
from pypulseq.opts import Opts
from pypulseq.calc_duration import calc_duration
from pypulseq.make_arbitrary_rf import make_arbitrary_rf
from pypulseq.make_sigpy_pulse import sigpy_n_seq
from pypulseq.sigpy_pulse_opts import SigpyPulseOpts
from dipy.io.gradients import read_bvals_bvecs

import spiraltraj
from sigpy.mri import rf as rfsig
import pulseq_helper as ph
from diffusion import diff_params
from prot import create_hdr
from verse import verse
from gre_refscan_B0 import gre_refscan_B0

#%% Parameters 
"""
PyPulseq units (SI): 
time:       [s] (not [ms] as in documentation)
spatial:    [m]
gradients:  [Hz/m] (gamma*T/m)
grad area:  [1/m]
flip angle: [rad]

Some units get converted below, others have to stay in non-SI units as spiral calculation needs different units.

SLR pulses should always be used for the diffusion sequence!!!
For multiband imaging, the use of VERSE pulses is recommended (minimum-time VERSE with energy constraint).
The pulse energy of the VERSE pulses will be printed. For a TR of 100ms, the combined energy of exc&ref should not exceed 100, if fatsat is used

Add new parameters also to the protocol parameter json file at the end of this script.
Use custom version of Pypulseq (important for delays): https://github.com/mavel101/pypulseq (branch dev_mv)
WIP: In Pulseq 1.4 delays are not used anymore (PyPulseq is still on 1.3.1 though)

Caution: Trapezoidal gradients will get a sign change in the physical coordinate system in old/compat rotation matrix mode

"""
# General
B0              = 2.893620   # field strength [T]
scanner         = 'connectom'   # scanner for acoustic resonance check - only 
grads_off       = False   # turn off gradients in reference scan (just for simulation of ECC)
seq_name        = 'axon_diameter' # sequence/protocol filename

# Sequence - Contrast and Geometry
fov             = 220         # field of view [mm]
TR              = 129.63         # repetition time [ms]
TE              = 51          # echo time [ms]
res             = 2.5         # in plane resolution [mm]
slice_res       = 2.5         # slice thickness [mm]
dist_fac        = 0      # distance factor for slices [%]
slices          = 54          # number of slices
averages        = 1           # number of averages
inner_avg       = True      # do averages in inner loop
repetitions     = 1         # number of repetitions

refscan         = 2       # 0: no refscan, 1: normal refscan, 2: B0 mapping refscan
res_refscan     = 1         # resolution of refscan, if B0 mapping is performed, 2mm is typically sufficient
bw_refscan      = 800       # Bandwidth of the reference scan [Hz]
flip_refscan    = 50        # reference scan flip angle
half_refscan    = False     # collect only half number of slices in refscan (with doubled slice thickness) - for 1mm datasets
separate_tr     = False     # separate TRs for echoes of B0 mapping refscans
prepscans       = 5          # number of preparation/dummy scans
noisescans      = 16          # number of noise scans

# ADC
os_factor       = 2         # oversampling factor (automatic 2x os from Siemens is not applied)
max_adc         = 16384       # maximum number of samples per ADC (originally 8192 was max, but VE12U seems to accept higher values)

# RF
flip_angle      = 90
slr_rf          = True       # True: use Sigpys SLR pulse design, False: use external or sinc pulses
verse_rf        = 2           # 0: no VERSE, 1: min SAR VERSE, 2: min time verse (recommended), 3: mintverse (hargreaves)
b1_max_fac      = [0.5,0.5]        # only verse_rf=2/3: factor to limit maximum B1 amplitude (in units of initial RF) [excitation, refocusing]
energy_fac      = [1,1]    # only verse_rf=2/3: factor to limit maximum energy (in units of initial RF) [excitation, refocusing]
rf_dur          = 5.2           # RF duration [ms]
tbp_exc         = 5           # time bandwidth product excitation pulse
rf_refoc_dur    = 9           # refocusing pulse duration [ms]
tbp_refoc       = 6           # time bandwidth product refocusing pulse
refoc_fac       = 1           # factor for slice thickness of refocusing pulse (=1 in Siemens diffusion seq)
rf_spoiling     = False       # RF spoiling
sms             = True        # multiband imaging?
sms_factor      = 2           # multiband factor
kz_steps        = 0           # number of steps in slice direction blipped spiral

fatsat          = True        # Fat saturation pulse
fatsat_tbp      = 2.1         # tbp of fatsat pulse [ms] (BW is fixed at 1000Hz atm) (increase at low TR to avoid SAR problems), DO NOT DECREASE, as pulse gets clipped

# Gradients
max_slew        = 182         # maximum slewrate [T/m/s] (system limit)
spiral_slew     = 120         # maximum slew rate of spiral gradients - for pre Emph: set lower than max_slew
max_grad        = 55          # maximum gradient amplitude [mT/m] (system limit)
max_grad_sp     = 42         # maximum gradient amplitude of spiral gradients - for pre_emph: set lower than max_grad

Nintl           = 2           # spiral interleaves
redfac          = 2           # reduction/acceleration factor
spiraltype      = 1           # 1: Spiral Out, 4: ROI, WIP: other spiral waveforms
spiral_os       = 1           # variable density spiral oversampling in center
trans_beg       = 0.33        # variable density transition beginning (between 0 and 1)
trans_end       = 0.4         # transition end (>trans_beg & between 0 and 1)

pre_emph        = False       # Gradient Pre-Emphasis based on GIRF
skope           = True       # add trigger for skope measurement
measure_delay   = False       # if False start skope measurement directly before spirals, if True start at the beginning of echo time delay
sync_scans      = 10          # number of Skope sync_scans   

# Diffusion
# number of b=0 images should be appr. directions/3.59 c.f. Kingsley: Intr. to DTI Part III (2006)
diff_slewrate   = 80          # diffusion gradient slewrate
diff_maxgrad    = 280         # diffusion gradient max grad strength
read_axon       = True        # read b-values and directions for axon diameter measurement (requires dipy)
b_val           = [0,700,1000] # b-values [s/mm^2] (b=0 has to be included)
directions      = [10,30,30]  # number of acquisitions for each b-value
delta           = 15          # duration of diffusion gradients [ms], spacing is automatically calculated
vol_TR          = None       # volume TR [s], if None take minimum volume TR

#%% Limits, checks and preparations

# Set System limits
rf_dead_time = 100e-6 # lead time before rf can be applied
rf_ringdown_time = 30e-6 # coil hold time (20e-6) + frequency reset time (10e-6)
system = Opts(max_grad=max_grad, grad_unit='mT/m', max_slew=max_slew, slew_unit='T/m/s', 
                rf_dead_time=rf_dead_time, rf_ringdown_time=rf_ringdown_time, grad_raster_time=ph.dt_grad, rf_raster_time=ph.dt_rf)

# convert parameters to Pulseq units
TR          *= 1e-3 # [s]
TE          *= 1e-3 # [s]
rf_dur      *= 1e-3 # [s]
rf_refoc_dur*= 1e-3 # [s]
slice_res   *= 1e-3 # [m]
res_refscan *= 1e-3 # [m]
delta       *= 1e-3 # [s]

# calculate effective interleaves
intl_eff = int(Nintl/redfac)

# averaging
if inner_avg:
    avgs_in = averages
    avgs_out = 1
else:
    avgs_in = 1
    avgs_out = averages

# do some checks
if sms:
    slices_eff = int(slices/sms_factor)
    slice_sep = slices/sms_factor*slice_res*(1+dist_fac*1e-2) # distance between multiband slices [m]
    if slices/sms_factor%1 != 0:
        raise ValueError('Number of slices is not multiple of sms factor')
    if slices/sms_factor%2 == 0:
        raise ValueError('Slices/sms_factor (= number of stacks) must be an odd number') # ref: Barth (08/2015)
    if sms_factor > 2:
        mb_phs = 'quad_mod' # saves SAR and peak amp
    else:
        mb_phs = 'None'
else:
    slices_eff = slices
    sms_factor = 1

if half_refscan and slices%2 != 0:
    raise ValueError('Reducing number of refscan slices only possible for even slice number.')

if skope:
    skope_delay = 200e-6 # delay/gradient free interval after Skope trigger
    min_dist = 200e-3 # minimum trigger distance due to relaxation
    trig_skip = int(np.ceil(min_dist/TR))
    if trig_skip > 1:
        print(f"TR too low to capture all triggers (minimum trigger distance 200ms). Only every {trig_skip}th trigger is captured.")
else:
    skope_delay = 0
    trig_skip = 0

if Nintl/redfac%1 != 0:
    raise ValueError('Number of interleaves is not multiple of reduction factor')

if int(fov/res+0.5) % 2:
    raise ValueError(f'Matrix size {int(fov/res+0.5)} (FOV/resolution) is not even.') 

if spiraltype!=1 and spiraltype!=4:
    ValueError('Right now only spiraltype 1 (spiral out) and 4 (ROI) possible.')

if redfac > 1 and refscan==0:
    print("WARNING: Cartesian reference scan is not activated.")

if len(b_val) != len(directions):
    raise ValueError('b-value and direction lists have to be the same length.')
if not b_val:
    raise ValueError('Select at least one b-value.')
if 0 not in b_val:
    print("WARNING: No b=0 axquisition selected.")
if rf_dur==rf_refoc_dur and tbp_exc==tbp_refoc:
    raise ValueError('Do not choose same duration and TBP for excitation and refocusing pulse. Crashes the sequence due to a Pulseq or Pypulseq bug.')

#%% RF Pulse and slab/slice selection gradient

# make rf pulse and calculate duration of excitation and rewinding
if slr_rf:
    sigpy_cfg = SigpyPulseOpts(pulse_type='slr', ptype='st')
    if flip_angle == 90:
        sigpy_cfg.ptype = 'ex'
        sigpy_cfg.cancel_alpha_phs = True
    rf, gz, gz_rew, rf_del = sigpy_n_seq(flip_angle=flip_angle*np.pi/180, system=system, duration=rf_dur, slice_thickness=slice_res,
                        time_bw_product=tbp_exc, pulse_cfg=sigpy_cfg, use='excitation', return_gz=True, return_delay = True, disp=False)
else:
    rf, gz, gz_rew, rf_del = make_sinc_pulse(flip_angle=flip_angle*np.pi/180, system=system, duration=rf_dur, slice_thickness=slice_res,
                            apodization=0.5, time_bw_product=tbp_exc, use='excitation', return_gz=True, return_delay = True)

if sms and not verse_rf:
    band_sep  = slice_sep/slice_res*tbp_exc # normalized distance between slices
    rf.signal = rfsig.multiband.mb_rf(rf.signal, n_bands=sms_factor, band_sep=band_sep, phs_0_pt=mb_phs)

if verse_rf:
    g_delay = 2e-6 # gradient delay, used for shifting RF pulse (approx. determined from Girf)
    vp_exc = verse(rf, gz, system, g_delay=g_delay)
    vp_exc.max_grad = min(2*gz.amplitude, system.max_grad) # limit the gradient amplitude - does this make sense?
    rf_max = np.max(abs(vp_exc.rf))
    rf_energy = vp_exc.calc_energy(use_input=True)
    print(f"RF energy before VERSE (excitation): {rf_energy:.2f}")
    if verse_rf == 1:
        vp_exc.verse_min_sar()
    elif verse_rf == 2:
        vp_exc.verse_min_time(B_max=b1_max_fac[0]*rf_max, E_max=energy_fac[0]*rf_energy)
    elif verse_rf == 3:
        vp_exc.mintverse(B_max=b1_max_fac[0]*rf_max, E_max=energy_fac[0]*rf_energy, ramp_grad=False)
    print(f"RF energy after VERSE (excitation): {vp_exc.calc_energy(use_input=False):.2f}")
    if sms:
        vp_exc.mb_verse(band_sep=slice_sep, n_bands=sms_factor, phs_0_pt=mb_phs)
    # the signal would be scaled incorrectly by "make_arbitrary_rf", therefore we calculate the signal below in the loop, where also the slice frequency modulation is applied
    rf, rf_del = make_arbitrary_rf(vp_exc.rf_vs, flip_angle=1, delay=vp_exc.ramp_len[0]*system.grad_raster_time, system=system, use='excitation', return_delay=True)
    rf_dur = rf.t[-1] # duration might have changed
    gz = make_arbitrary_grad(channel='z', waveform=vp_exc.g_vs, system=system)
    if rf.delay > vp_exc.ramp_len[0]*system.grad_raster_time:
        gz.delay = rf.delay - vp_exc.ramp_len[0]*system.grad_raster_time
    gz_rew = make_trapezoid(channel='z', area=-1*np.sum(gz.waveform)/2*system.grad_raster_time, system=system)        

# refocusing pulse
if slr_rf:
    sigpy_cfg_ref = SigpyPulseOpts(pulse_type='slr', ptype='se')
    rf_refoc, gz_refoc, _ = sigpy_n_seq(flip_angle=2*flip_angle*np.pi/180, system=system, duration=rf_refoc_dur, slice_thickness=slice_res*refoc_fac,
                        time_bw_product=tbp_refoc, pulse_cfg=sigpy_cfg_ref, use='refocusing', return_gz=True, disp=False)
else:
    rf_refoc, gz_refoc, _  = make_sinc_pulse(flip_angle=2*flip_angle*np.pi/180, system=system, duration=rf_refoc_dur, slice_thickness=slice_res*refoc_fac,
                            apodization=0.5, time_bw_product=tbp_refoc, use='refocusing', return_gz=True)
gz_refoc_area = gz_refoc.area

if sms and not verse_rf:
    band_sep_refoc  = slice_sep/slice_res*tbp_refoc
    rf_refoc.signal = rfsig.multiband.mb_rf(rf_refoc.signal, n_bands=sms_factor, band_sep=band_sep_refoc, phs_0_pt=mb_phs)

if verse_rf:
    vp_ref = verse(rf_refoc, gz_refoc, system, g_delay=g_delay)
    vp_ref.max_grad = min(2*gz_refoc.amplitude, system.max_grad)
    rf_ref_max = np.max(abs(vp_ref.rf))
    rf_ref_energy = vp_ref.calc_energy(use_input=True)
    print(f"RF energy before VERSE (refocusing): {rf_ref_energy:.2f}")
    if verse_rf == 1:
        vp_ref.verse_min_sar()
    elif verse_rf == 2:
        vp_ref.verse_min_time(B_max=b1_max_fac[1]*rf_ref_max, E_max=energy_fac[1]*rf_ref_energy)
    elif verse_rf == 3:
        vp_ref.mintverse(B_max=b1_max_fac[1]*rf_ref_max, E_max=energy_fac[1]*rf_ref_energy, ramp_grad=False)
    print(f"RF energy after VERSE (refocusing): {vp_ref.calc_energy(use_input=False):.2f}")
    if sms:
        vp_ref.mb_verse(band_sep=slice_sep, n_bands=sms_factor, phs_0_pt=mb_phs)
    rf_refoc = make_arbitrary_rf(vp_ref.rf_vs, flip_angle=1, system=system, use='refocusing')
    gz_refoc = make_arbitrary_grad(channel='z', waveform=vp_ref.g_vs, system=system)
    gz_refoc_area = np.sum(gz_refoc.waveform)*system.grad_raster_time

# define crusher - only for b=0
#  rf pulse                  #########
#  slice gradient          #############
#  crusher gradients   ######         ######

system.max_slew = 90 * system.gamma # avoid stimulation
system.max_grad = 45e-3 * system.gamma

crusher_area = gz_refoc_area / 2 # define as spoiler
amp_crush, ftop_crush, ramp_crush = ph.trap_from_area(crusher_area, system)
crusher_z1 = make_trapezoid(channel='z',system=system, amplitude=amp_crush, flat_time=ftop_crush, rise_time=ramp_crush)
crusher_z2 = copy(crusher_z1)
crusher_dur = calc_duration(crusher_z1)

system.max_slew = max_slew * system.gamma
system.max_grad = 1e-3*max_grad * system.gamma

# merge crushers with refocusing gradient
rf_refoc.delay = crusher_dur
if verse_rf:
    grad_refoc = ph.merge_verse_ramps(vp_ref, crusher_z1, system)
    grad_refoc = make_arbitrary_grad(channel='z', waveform=grad_refoc, system=system)
else:
    gz_refoc = make_trapezoid(channel='z', system=system, amplitude=gz_refoc.amplitude, flat_time=gz_refoc.flat_time, rise_time=ramp_crush)
    grad_refoc = ph.merge_ramps([crusher_z1, gz_refoc, crusher_z2], system=system) # cumulative area is preserved
refoc_dur  = calc_duration(grad_refoc)

# for b=0 define crushers also on x- and y channels
crusher_xy1 = copy(crusher_z1)
crusher_xy2 = copy(crusher_z1)
crusher_xy2.delay = round(calc_duration(grad_refoc) - crusher_dur, ndigits=5)
crusher_xy1.channel = 'x'
crusher_xy2.channel = 'x'
crusher_x = ph.add_gradients([crusher_xy1,crusher_xy2], system=system)
crusher_xy1.channel = 'y'
crusher_xy2.channel = 'y'
crusher_y = ph.add_gradients([crusher_xy1,crusher_xy2], system=system)

# timing variables
exc_to_rew = calc_duration(rf, gz, rf_del) - rf.delay - ph.round_up_to_raster(rf_dur/2, decimals=5) # time from middle of rf pulse to rewinder
rew_dur = calc_duration(gz_rew)

# RF spoiling parameters
rf_spoiling_inc = 50 # increment of RF spoiling [°]
rf_phase        = 0 
rf_inc          = 0

# Fat saturation
if fatsat:
    fatsat_fa = 110 # flip angle [°]
    fw_shift = 3.35e-6 # unsigned fat water shift [ppm]
    offset = -1 * int(B0*system.gamma*fw_shift)
    fatsat_bw = abs(offset) # bandwidth [Hz]
    fatsat_dur = ph.round_up_to_raster(fatsat_tbp/fatsat_bw, decimals=5)
    rf_fatsat, fatsat_del = make_gauss_pulse(flip_angle=fatsat_fa*np.pi/180, duration=fatsat_dur, bandwidth=fatsat_bw, freq_offset=offset, system=system, return_delay = True)

#%% Diffusion gradients and delay calculation 

# Delays in diffusion sequence

#  RF  rew  diffgrad1   refoc_delay   RF_refoc   spacing delay  diffgrad2   te_delay   readout
# #### ###  ##########     ######       ####       ######      ##########   ######   ############
#   -------------- TE/2 ------------------
#           ----------------    spacing   ---------------------
#   -------------------------------------------- TE ---------------------------------


diff_maxgrad_Hz = 1e-3 * diff_maxgrad * system.gamma
diff_slewrate_Hz = diff_slewrate * system.gamma
system.max_slew = diff_slewrate_Hz
system.max_grad = diff_maxgrad_Hz

# read b-values & directions
axon_file = "axon_diameter"
bval_list, dir_list = read_bvals_bvecs(axon_file + ".bval", axon_file + ".bvec")
n_dirs = len(dir_list)

# calculate diffusion parameters
b_val_max = max(bval_list)
params_diff = diff_params(b_val=b_val_max, delta=delta, spacing=None)
params_diff.calc_params(system)

diff_list = []
for k, bval in enumerate(bval_list):
    diff_list.append({"bval": int(bval), "dir": dir_list[k]})

# only for saving
b_val = list(set(bval_list))
directions = [np.count_nonzero(bval_list==b_val[k]) for k in range(len(b_val))]

# calculate diffusion gradient duration
diffgrad_ramp = ph.round_up_to_raster(system.max_grad/system.max_slew, 5)
diffgrad_flat = params_diff.delta - diffgrad_ramp
diffgrad_dur = diffgrad_flat + 2*diffgrad_ramp

if params_diff.spacing < refoc_dur+diffgrad_dur:
    raise ValueError('Spacing between diffusion gradients too low. Increase spacing or decrease duration (delta) of diffusion gradients.')

# make diffusion gradients - x-channel will get sign change from Pulseq rotation matrix (in the old/compat Pulseq mode), if gradient is trapezoidal
trap_diff_x = make_trapezoid(channel='x', system=system, flat_time=diffgrad_flat, rise_time=diffgrad_ramp, amplitude=-1*system.max_grad, max_slew=system.max_slew)
trap_diff_y = make_trapezoid(channel='y', system=system, flat_time=diffgrad_flat, rise_time=diffgrad_ramp, amplitude=system.max_grad, max_slew=system.max_slew)
trap_diff_z = make_trapezoid(channel='z', system=system, flat_time=diffgrad_flat, rise_time=diffgrad_ramp, amplitude=system.max_grad, max_slew=system.max_slew)
diff_gradients = [trap_diff_x,trap_diff_y,trap_diff_z]

# reset the maximum slewrate and gradient
system.max_slew = max_slew * system.gamma
system.max_grad = 1e-3 * max_grad * system.gamma

# calculate minimum and maximum TE
te1 = 2*(exc_to_rew + rew_dur + diffgrad_dur + ph.round_up_to_raster(refoc_dur/2, decimals=5)) # min TE if spacing is small
te2 = exc_to_rew + rew_dur + params_diff.spacing + diffgrad_dur # min TE if spacing is large
min_te = max(te1,te2) # min TE: readout follows directly after 2nd diffusion gradient
if skope:
    min_te += skope_delay
if TE < min_te:
    raise ValueError(f'Minimum TE: {min_te*1e3}')

# calculate delay after slice rewinder
# if the refocusing pulse is directly before 2nd diff gradient, we have to add an additional delay to make the TE possible
max_te = 2*(exc_to_rew + rew_dur + params_diff.spacing - ph.round_up_to_raster(refoc_dur/2, decimals=5))
if max_te < TE:
    extra_delay = TE-max_te
    diff_delay = make_delay(d=extra_delay/2)
else:
    diff_delay = make_delay(d=0)

# make delays
refoc_delay = make_delay(d = TE/2 - diff_delay.delay - exc_to_rew - rew_dur - diffgrad_dur - ph.round_up_to_raster(refoc_dur/2, decimals=5))
spacing_delay = make_delay(d = params_diff.spacing - diffgrad_dur - refoc_delay.delay - refoc_dur)
te_delay = make_delay(d = TE/2 - ph.round_up_to_raster(refoc_dur/2, decimals=5) - spacing_delay.delay - diffgrad_dur) # delay to achieve protocol TE


#%% Spiral Readout Gradients

# Parameters spiral trajectory:

# parameter         description               default value
# ---------        -------------              --------------

# nitlv:      number of spiral interleaves        15
# res:        resolution                          1 mm
# fov:        target field of view                192 mm
# max_amp:    maximum gradient amplitude          42 mT/m
# min_rise:   minimum gradient risetime           5 us/(mT/m)
# spiraltype: 1: spiral out                   
#             2: spiral in                        
#             3: double spiral                    x
#             4: ROI
#             5: RIO
# spiral_os:  spiral oversampling in center       1

# Maximum rotation angle for spirals
if spiraltype==3:
    max_rot     = np.pi
else:
    max_rot     = 2*np.pi  

# read in Spirals [T/m]
min_rise_sp = 1/spiral_slew * 1e3
spiral_calc = spiraltraj.calc_traj(nitlv=Nintl, fov=fov, res=res, spiraltype=spiraltype,
                             min_rise=min_rise_sp, max_amp=max_grad_sp, spiral_os=spiral_os, 
                             vd_transition_begin=trans_beg, vd_transition_end=trans_end)
spiral_calc = np.asarray(spiral_calc)
spiral_x = 1e-3*spiral_calc[:,0]
spiral_y = 1e-3*spiral_calc[:,1]

N_spiral = len(spiral_x)
readout_dur = N_spiral*system.grad_raster_time # readout duration [s]

# write spiral readout blocks to list
spirals = [{'deph': [None, None], 'spiral': [None, None], 'reph': [None, None]} for k in range(Nintl)]
reph_dur = []
save_sp = np.zeros((Nintl, 2, N_spiral)) # save gradients for FIRE reco
rot_angle = np.linspace(0, max_rot, Nintl, endpoint=False)
for k in range(Nintl):
    # rotate spiral gradients for shot selection
    sp_x, sp_y = ph.rot_grad(spiral_x, spiral_y, rot_angle[k])

    save_sp[k,0,:] = sp_x
    save_sp[k,1,:] = sp_y

    # unit to [Hz/m], make spiral gradients
    sp_x *= system.gamma
    sp_y *= system.gamma

    spiral_delay = 20e-6 # delay to avoid ADC artifact (first few points of ADC might be corrupted)
    spirals[k]['spiral'][0] = make_arbitrary_grad(channel='x', waveform=sp_x, delay=spiral_delay, system=system)
    spirals[k]['spiral'][1] = make_arbitrary_grad(channel='y', waveform=sp_y, delay=spiral_delay, system=system)

    if spiraltype==1:
        # calculate rephaser area
        area_x = sp_x.sum()*system.grad_raster_time
        area_y = sp_y.sum()*system.grad_raster_time

        # calculate rephasers and make gradients - add spoiler area to rephaser
        spoiler_area = np.sum(ph.waveform_from_seqblock(gz))*system.grad_raster_time/np.sqrt(3)
        amp_x, ftop_x, ramp_x = ph.trap_from_area(-area_x+spoiler_area, system, slewrate=100, max_grad=30e-3) # reduce slew rate & max_grad to to avoid stimulation
        amp_y, ftop_y, ramp_y = ph.trap_from_area(-area_y+spoiler_area, system, slewrate=100, max_grad=30e-3)
        spirals[k]['reph'][0] = make_trapezoid(channel='x', system=system, amplitude=amp_x, flat_time=ftop_x, rise_time=ramp_x)
        spirals[k]['reph'][1] = make_trapezoid(channel='y', system=system, amplitude=amp_y, flat_time=ftop_y, rise_time=ramp_y)
        reph_dur.append(max(ftop_x+2*ramp_x, ftop_y+2*ramp_y))

# check spirals for acoustic resonances
if B0 > 4:
    resonances = [(500,600), (930, 1280)] # 7T resonances
else:
    if scanner == 'skyra':
        resonances = [(535, 635), (1010,1230)] # 3T Skyra resonances
    elif scanner == 'connectom':
        resonances = [(280,340), (546, 646), (1000,1500)] # 3T Connectom resonances
    else:
        raise ValueError('Unknown scanner name for 3T, select either skyra or connectom.')
freq_max = ph.check_resonances([spiral_x,spiral_y], resonances) 

#%% Spoiler Gradients

spoiler_area = np.sum(ph.waveform_from_seqblock(gz))*system.grad_raster_time/np.sqrt(3) # moment of excitation pulse
amp_spoil, ftop_spoil, ramp_spoil = ph.trap_from_area(spoiler_area, system, slewrate=100, max_grad = 30e-3) # reduce slew rate and max_grad to avoid stimulation

spoiler_x = make_trapezoid(channel='x',system=system, amplitude=amp_spoil, flat_time=ftop_spoil, rise_time=ramp_spoil)
spoiler_y = make_trapezoid(channel='y',system=system, amplitude=amp_spoil, flat_time=ftop_spoil, rise_time=ramp_spoil)
spoiler_z = make_trapezoid(channel='z',system=system, amplitude=amp_spoil, flat_time=ftop_spoil, rise_time=ramp_spoil)
spoiler_x_neg  = make_trapezoid(channel='x',system=system, amplitude=-1*amp_spoil, flat_time=ftop_spoil, rise_time=ramp_spoil)
spoiler_y_neg  = make_trapezoid(channel='y',system=system, amplitude=-1*amp_spoil, flat_time=ftop_spoil, rise_time=ramp_spoil)
spoiler_z_neg  = make_trapezoid(channel='z',system=system, amplitude=-1*amp_spoil, flat_time=ftop_spoil, rise_time=ramp_spoil)

spoiler_dur = calc_duration(spoiler_z)

#%% ADC

max_grad_sp_cmb = 1e3*np.max(np.sqrt(abs(spiral_x)**2+abs(spiral_y)**2))
dwelltime = 1/(system.gamma*max_grad_sp_cmb*fov*os_factor)*1e6 # ADC dwelltime [s]
dwelltime = ph.trunc_to_raster(dwelltime, decimals=7) # truncate dwelltime to 100 ns (scanner limit)
min_dwelltime = 1e-6
if dwelltime < min_dwelltime:
    dwelltime = min_dwelltime
print(f"ADC dwelltime: {1e6*dwelltime} us")

num_samples = round((readout_dur+spiral_delay)/dwelltime)
if num_samples%2==1:
    num_samples += 1 # even number of samples

if num_samples <= max_adc:
    num_segments = 1
    print('Number of ADCs: {}.'.format(num_samples))
else:
    # the segment duration has to be on the gradient raster
    # increase number of segments or samples/segments to achieve this
    # number of samples and number of samples per segment should always be an even number
    num_segments = 2
    if (num_samples/num_segments % 2 != 0):
        num_samples += 2
    segm_dur = 1e5 * dwelltime * num_samples/num_segments # segment duration [10us - gradient raster]
    while (not round(segm_dur,ndigits=5).is_integer() or num_samples/num_segments > 8192):
        if num_samples/num_segments > 8192:
            num_segments += 1
            while (num_samples/num_segments % 2 != 0):
                num_samples += 2
        else:
            num_samples += 2*num_segments
        segm_dur = 1e5 * dwelltime * num_samples/num_segments 
    print('ADC has to be segmented!! Number of ADCs: {}. Per segment: {}. Segments: {}.'.format(num_samples,num_samples/num_segments,num_segments))

    # self check
    if (num_samples/num_segments % 2 != 0 or num_samples % 2 != 0 or not round(segm_dur,ndigits=5).is_integer()):
        raise ValueError("Check if number of samples and number of samples per segment are even. Check if segment duration is on gradient raster time.")

if num_samples > 65535: # max of uint16 used by ISMRMRD
    raise ValueError("Too many samples for ISMRMRD format - lower the oversampling factor or take more interleaves")

adc = make_adc(system=system, num_samples=num_samples, dwell=dwelltime)
adc_dur = calc_duration(adc)
adc_delay = ph.round_up_to_raster(adc_dur+200e-6, decimals=5) # add small delay after readout for ADC frequency reset event and to avoid stimulation by rephaser
adc_delay = make_delay(d=adc_delay)
if skope:
    t_skope = (adc_dur+te_delay.delay+1e-3)*1e3 # add 1 ms to be safe
    print('Minimum Skope acquisition time: {:.2f} ms'.format(t_skope))

#%% Set up protocol for FIRE reco and write header

date = datetime.date.today().strftime('%Y%m%d')
filename = date + '_' + seq_name

# set some parameters for the protocol
t_min = dwelltime/2

# create new directory if needed
ismrmrd_file = f"{filename}.h5"

# set up protocol file and create header
if os.path.exists(ismrmrd_file):
    raise ValueError("Protocol name already exists. Choose different name")
prot = ismrmrd.Dataset(ismrmrd_file)
hdr = ismrmrd.xsd.ismrmrdHeader()
params_hdr = {"trajtype": "spiral", "fov": fov, "fov_z": slice_res, "res": res, "slices": slices, "slice_res": slice_res*(1+dist_fac*1e-2), 
                "nintl": intl_eff, "avg": averages, "rep": repetitions, "ncontrast": n_dirs,
                "nsegments": num_segments, "dwelltime": dwelltime, "traj_delay": spiral_delay, "t_min": t_min, 
                "os_region": trans_beg, "os_factor": os_factor, "redfac": redfac, "sms_factor": sms_factor, "half_refscan": half_refscan}
create_hdr(hdr, params_hdr)

#%% Add sequence blocks to sequence & write acquisitions to protocol

# Set up the sequence
seq = Sequence()
trig_ctr = 0

# Definitions section in seq file
seq.set_definition("Name", filename) # protocol name is saved in Siemens header for FIRE reco
seq.set_definition("FOV", [1e-3*fov, 1e-3*fov, slice_res*(1+dist_fac*1e-2)*(slices-1)+slice_res]) # this sets the volume display in the UI
seq.set_definition("Slice_Thickness", "%f" % slice_res) # this sets the receive gain
if num_segments > 1:
    seq.set_definition("MaxAdcSegmentLength", "%d" % int(num_samples/num_segments+0.5)) # for automatic ADC segment length setting

# TokTokTok
tokx = make_trapezoid(channel='x', amplitude=-1e-3*system.gamma, rise_time=1e-3, flat_time=4e-3)
toky = make_trapezoid(channel='y', amplitude=1e-3*system.gamma, rise_time=1e-3, flat_time=4e-3)
tokz = make_trapezoid(channel='z', amplitude=1e-3*system.gamma, rise_time=1e-3, flat_time=4e-3)
seq.add_block(tokx,toky,tokz,make_delay(d=0.5))
seq.add_block(tokx,toky,tokz,make_delay(d=0.5))
seq.add_block(tokx,toky,tokz,make_delay(d=0.5))

# Noise scans
noise_samples = 256
noise_adc = make_adc(system=system, num_samples=256, dwell=dwelltime, delay=10e-6) # delay to be safe with pTx system (had crashes due to short NCO gaps)
noise_delay = make_delay(d=ph.round_up_to_raster(calc_duration(noise_adc)+1e-3,decimals=5)) # add some more time to the ADC delay to be safe
for k in range(noisescans):
    seq.add_block(noise_adc, noise_delay)
    acq = ismrmrd.Acquisition()
    acq.setFlag(ismrmrd.ACQ_IS_NOISE_MEASUREMENT)
    prot.append_acquisition(acq)

# Perform cartesian reference scan
if refscan:
    te_1 = 2.42e-3
    te_2 = 4.84e-3
    te_3 = 4.84e-3
    rf_ref_dur = 1.2e-3
    tbp_ref_dur = 4
    if refscan == 1:
        ref_lines = 30 # ecalib takes 24 lines as default
        center_out = False # center out readout
        TE_ref = [te_1]
        params_ref = {"TE": TE_ref, "fov":fov*1e-3, "ref_lines": ref_lines, "slices":slices, "slice_res":slice_res, "dist_fac": dist_fac,
                    "flip_angle":flip_refscan, "readout_bw": bw_refscan, 'center_out': center_out, "rf_dur": rf_ref_dur, "tbp": tbp_ref_dur}
        gre_refscan_B0(seq, prot=prot, system=system, params=params_ref, grads_off=grads_off)
    elif refscan == 2:
        if half_refscan:
            slices_ref = slices // 2
            slice_res_ref = slice_res * 2
        else:
            slices_ref = slices
            slice_res_ref = slice_res
        center_out = True # center out readout
        if separate_tr:
            TE_ref = [te_1, te_2]
        else:
            TE_ref = [te_1, te_3]
        params_ref = {"TE": TE_ref, "fov":fov*1e-3, "res":res_refscan, "slices":slices_ref, "slice_res":slice_res_ref, "dist_fac": dist_fac,
                    "flip_angle":flip_refscan, "readout_bw": bw_refscan, 'separate_tr': separate_tr, 'center_out': center_out, "rf_dur": rf_ref_dur, "tbp": tbp_ref_dur}
        gre_refscan_B0(seq, prot=prot, system=system, params=params_ref, grads_off=grads_off)
    else:
        raise ValueError("Invalid refscan selection.")

dur_until_ref = seq.duration()[0]
print(f"Sequence duration after reference scan: {dur_until_ref:.2f} s")

# Skope sync scans
if skope:
    if measure_delay:
        trig_delay = 0 # measure at the beginning of the echo time delay
    else:
        trig_delay = te_delay.delay - skope_delay # measure 200us before spiral readout

    adc_sync = make_adc(system=system, num_samples=4000, dwell=dwelltime)
    adc_sync_delay = make_delay(d=ph.round_up_to_raster(calc_duration(adc_sync)+200e-6, decimals=5))
    trig = make_digital_output_pulse(channel='ext1', duration=system.grad_raster_time, delay=trig_delay)

    for j in range(sync_scans):
        seq.add_block(trig, te_delay)
        seq.add_block(adc_sync, adc_sync_delay)
        seq.add_block(make_delay(d=50e-3)) # some delay between triggers

        acq = ismrmrd.Acquisition()
        acq.setFlag(ismrmrd.ACQ_IS_DUMMYSCAN_DATA)
        prot.append_acquisition(acq)

    sync_scan_delay = make_delay(d=5)
    seq.add_block(sync_scan_delay) # Skope trigger receipt has dead time after sync scans

else:
    sync_scans = 0 

""" diffusion

The following code generates a Pulseq diffusion sequence.
Single-shot & multishot acquisitions are possible
For multishot, sufficient oversampling in kspace center has to be chosen, as a phase correction is needed.

"""

vol_TR_delay = vol_TR - (TR*slices_eff) if vol_TR is not None else None
if vol_TR is None:
    vol_TR = TR*slices_eff
print(f"Volume TR: {vol_TR:.3f} s.")

# save loop variables before prepscans
repetitions_ = repetitions
avgs_in_ = avgs_in
avgs_out_ = avgs_out
intl_eff_ = intl_eff
diff_list_ = diff_list.copy()

# run prepscans, then imaging scans
for prep in range(prepscans+1):
    if prep != prepscans:
        repetitions = avgs_out = avgs_in = intl_eff = 1
        diff_list = [{"bval": 0, "dir": np.zeros(3)}]
    else:
        repetitions = repetitions_
        avgs_in = avgs_in_
        avgs_out = avgs_out_
        intl_eff = intl_eff_
        diff_list = diff_list_.copy()

    for rep in range(repetitions):
        for avg_out in range(avgs_out):
            for vol_ix, diff in enumerate(diff_list):
                diff_gradients[0].amplitude = -1 * np.sqrt(diff['bval']/b_val_max) * diff['dir'][0] * diff_maxgrad_Hz # consider sign change from Pulseq rotation matrix
                diff_gradients[1].amplitude = np.sqrt(diff['bval']/b_val_max) * diff['dir'][1] * diff_maxgrad_Hz
                diff_gradients[2].amplitude = np.sqrt(diff['bval']/b_val_max) * diff['dir'][2] * diff_maxgrad_Hz
                for avg_in in range(avgs_in):
                    for n in range(intl_eff):
                        # slice ordering acc to Siemens method
                        if slices_eff%2 == 1:
                            slc = 0
                        else:
                            slc = 1
                        for slc_ctr in range(slices_eff):
                            if slc_ctr==int(slices_eff/2+0.5):
                                if slices_eff%2 == 1:
                                    slc = 1
                                else:
                                    slc = 0

                            # Add (fatsat and) excitation pulse
                            if fatsat:
                                rf_fatsat.phase_offset = rf_phase / 180 * np.pi # always use RF spoiling for fat sat pulse
                                seq.add_block(rf_fatsat, fatsat_del)
                                rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
                                rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]
                                seq.add_block(spoiler_x_neg, spoiler_y_neg, spoiler_z_neg)

                            if verse_rf:
                                rf.signal = vp_exc.rf_vs * np.exp(1j* 2*np.pi * vp_exc.mod * (slice_res * (slc - (slices_eff - 1) / 2) * (1+dist_fac*1e-2)))
                                rf_refoc.signal = vp_ref.rf_vs * np.exp(1j* 2*np.pi * vp_ref.mod * (slice_res * (slc - (slices_eff - 1) / 2) * (1+dist_fac*1e-2)))
                            else:
                                rf.freq_offset = gz.amplitude * slice_res * (slc - (slices_eff - 1) / 2) * (1+dist_fac*1e-2)
                                rf_refoc.freq_offset = gz_refoc.amplitude * slice_res * (slc - (slices_eff - 1) / 2) * (1+dist_fac*1e-2)
                            seq.add_block(rf,gz,rf_del)
                            seq.add_block(gz_rew)

                            # diffusion block
                            seq.add_block(diff_delay)
                            seq.add_block(diff_gradients[0],diff_gradients[1],diff_gradients[2])
                            seq.add_block(refoc_delay)
                            if diff['bval'] == 0:
                                seq.add_block(rf_refoc,grad_refoc,crusher_x,crusher_y)
                            else:
                                seq.add_block(rf_refoc,grad_refoc)
                            seq.add_block(spacing_delay)
                            seq.add_block(diff_gradients[0],diff_gradients[1],diff_gradients[2])

                            # Skope trigger - keep minimum distance of 200us between subsequent triggers
                            if skope and slc_ctr%trig_skip==0 and slices_eff-slc_ctr >= trig_skip and prep == prepscans:
                                trig_ctr += 1
                                seq.add_block(trig, te_delay)
                            else:
                                seq.add_block(te_delay)

                            # spiral readout block
                            spiral_block = [spirals[n*redfac]['spiral'][0], spirals[n*redfac]['spiral'][1], adc_delay]
                            if prep == prepscans:
                                spiral_block.append(adc)
                            seq.add_block(*spiral_block)

                            # delay for TR
                            if spiraltype==1:
                                min_tr = rf.delay + ph.round_up_to_raster(rf_dur/2, decimals=5) + TE + adc_delay.delay + calc_duration(spirals[n*redfac]['reph'][0],spirals[n*redfac]['reph'][1], spoiler_z)
                            else:
                                min_tr = rf.delay + ph.round_up_to_raster(rf_dur/2, decimals=5) + TE + adc_delay.delay + spoiler_dur
                            if fatsat:
                                min_tr += fatsat_del.delay + spoiler_dur
                            if TR < min_tr:
                                raise ValueError('Minimum TR is {} ms.'.format(min_tr*1e3))
                            tr_delay = make_delay(d=TR-min_tr)
                            seq.add_block(tr_delay)

                            # spoiler/rephaser directly before next fatpulse
                            if spiraltype==1:
                                seq.add_block(spirals[n*redfac]['reph'][0],spirals[n*redfac]['reph'][1], spoiler_z)
                            else:
                                seq.add_block(spoiler_x, spoiler_y, spoiler_z)

                            # add protocol information
                            if prep == prepscans:
                                for seg in range(num_segments):
                                    acq = ismrmrd.Acquisition()
                                    if (n == intl_eff-1) and (seg == num_segments-1):
                                        acq.setFlag(ismrmrd.ACQ_LAST_IN_SLICE)
                                    acq.idx.kspace_encode_step_1 = n
                                    acq.idx.slice = slc
                                    acq.idx.contrast = vol_ix
                                    acq.idx.average = max(avg_in, avg_out)
                                    acq.idx.repetition = rep
                                    acq.idx.segment = seg
                                    acq.user_int[0] = diff['bval']
                                    acq.user_float[:3] = diff['dir']
                                    
                                    # save gradient only in first segment to save space
                                    if seg == 0:
                                        # use the trajectory field for the gradient array
                                        acq.resize(trajectory_dimensions = save_sp.shape[1], number_of_samples=save_sp.shape[2], active_channels=0)
                                        acq.traj[:,:2] = np.swapaxes(save_sp[n*redfac],0,1) # [samples, dims]
                                    prot.append_acquisition(acq)
                            
                            slc += 2 # interleaved slice acquisition

                if vol_TR_delay is not None:
                    seq.add_block(make_delay(d=vol_TR_delay))

                        # slices
                    # intl
                # avg_in
            # contrast
        # avg_out
    # reps
# prepscans

# save b-values and directions as array
prot.append_array("b_values", np.asarray(bval_list, dtype=np.float32))
prot.append_array("Directions", np.asarray(dir_list, dtype=np.float32))

if skope:
    print(f"Number of Skope triggers: {trig_ctr}.")

#%% write sequence & protocol

seq.write(f'{filename}.seq')
prot.close()
