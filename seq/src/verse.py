"""
VERSE pulse class

"""

import numpy as np
import types
import pypulseq as pp

def ifft(sig, dim=None):
    """ Computes the Fourier transform from k-space to image space 
    along a given or all dimensions

    :param img: image space data
    :param dim: vector of dimensions to transform
    :returns: data in k-space (along transformed dimensions)
    """
    import collections.abc
    if dim is None:
        dim = range(sig.ndim)
    elif not isinstance(dim, collections.abc.Iterable):
        dim = [dim]

    sig = np.fft.ifftshift(sig, axes=dim)
    sig = np.fft.ifftn(sig, axes=dim)
    sig = np.fft.fftshift(sig, axes=dim)

    return sig

def fft(sig, dim=None):
    """ Computes the Fourier transform from image space to k-space
    along a given or all dimensions

    :param img: image space data
    :param dim: vector of dimensions to transform
    :returns: data in k-space (along transformed dimensions)
    """
    import collections.abc
    if dim is None:
        dim = range(sig.ndim)
    elif not isinstance(dim, collections.abc.Iterable):
        dim = [dim]

    sig = np.fft.ifftshift(sig, axes=dim)
    sig = np.fft.fftn(sig, axes=dim)
    sig = np.fft.fftshift(sig, axes=dim)

    return sig

class verse():

    def __init__(self, rf, g, system, g_delay=0.0):
        """
        VERSE pulse class

        input:
        rf: Base non-VERSE RF waveform (Pypulseq object or Array) on RF raster time (typically 1us)
        gz: corresponding slice gradient (Pypulseq object or value with gradient amplitude [Hz] or numpy array on same raster as RF)
        system: PyPulseq system object or dict containing system parameters (max_grad, max_slew, rf_raster_time, grad_raster_time)
        g_delay: gradient delay [s] - delays are are around 2us-3us for x/y and around 0us for z

        units are [Hz] and [Hz/m] as in PyPulseq

        Main functions are verse_min_sar and verse_min time, which calculate SAR or time-optimized VERSE pulses

        Outputs rf_vs and g_vs are on RF and gradient raster time respectively
        """

        # Input
        if type(rf) ==  types.SimpleNamespace:
            self.rf = rf.signal
        else:
            self.rf = rf
        self.N = len(self.rf)

        if type(g) ==  types.SimpleNamespace:
            self.G = np.ones(self.N) * g.amplitude
        elif isinstance(g,np.ndarray):
            self.G = g
        else:
            self.G = np.ones(self.N) * g
        self.g_delay = g_delay

        if any(self.G<0):
            raise ValueError("Gradient values have to be non-negative.") 

        if type(system) ==  pp.opts.Opts:
            self.max_grad = system.max_grad
            self.max_slew = system.max_slew
            self.grad_raster_time = system.grad_raster_time
            self.rf_raster_time = system.rf_raster_time
        else:
            self.max_grad = system['max_grad']
            self.max_slew = system['max_slew']
            self.grad_raster_time = system['grad_raster_time']
            self.rf_raster_time = system['rf_raster_time']

        # Output
        self.rf_vs = None
        self.g_vs = None
        self.N_vs = self.N
        self.mod = None # modulation for slice offset [1/m]
        self.ramp_len = (0,0) # ramp length of g_vs

        # Helper
        self._t_vs = None
        self._alpha = None 
        self.safety_fac = 0.99 # safety factor to avoid slew rate overflow
        self.max_slew *= self.safety_fac

    def verse_min_sar(self, add_ramps=True, use_girf=False):
        """
        Calculate minimum-SAR VERSE pulses from input rf pulse (Conolly 1988)
        Calculations are done on rf raster time (1us) and waveforms are interpolated in the end

        use_girf: Use girf for modulation function calculation
        """

        # calculate minimum-SAR VERSE gradient
        c = self.G/self.N * np.sum(abs(self.rf))
        g = np.divide(c, abs(self.rf), out=np.ones_like(self.rf)*self.max_grad, where=abs(self.rf)!=0)
        ix_h = np.nonzero(g >= self.max_grad)
        n_h = len(ix_h[0])
        g[ix_h] = self.max_grad
        ix_l = np.nonzero(g < self.max_grad)
        not_cvgd = True
        while not_cvgd:
            n_h_bf = n_h
            c = np.sum(abs(self.rf[ix_l])) / (self.N/self.G - n_h/self.max_grad)
            g[ix_l] = c / abs(self.rf[ix_l])
            ix_h = np.nonzero(g >= self.max_grad)
            n_h = len(ix_h[0])
            g[ix_h] = self.max_grad
            ix_l = np.nonzero(g < self.max_grad)
            not_cvgd = bool(n_h-n_h_bf)

        # Calculate waveforms 
        self.g_vs = g
        self._alpha = self.g_vs / self.G
        self.rf_vs = self.rf * self._alpha
        self._t_vs = self.rf_raster_time/self._alpha

        # adress possible slewrate overflow
        self._smooth_gradient()

        # interpolate to rf and gradient raster
        raster_ratio = round(self.grad_raster_time/self.rf_raster_time)
        vs_raster = np.cumsum(self._t_vs)
        rf_raster =  np.arange(0.5,self.N_vs+0.5) * self.rf_raster_time
        g_raster = np.arange(0.5,self.N_vs//raster_ratio+0.5) * self.grad_raster_time
        self.rf_vs = self._intp_wf(self.rf_vs, vs_raster, rf_raster)
        self.g_vs = self._intp_wf(self.g_vs, vs_raster, g_raster)
        
        print(f"VERSE pulse length: {1e3*self.N_vs*self.rf_raster_time} ms")

        # delay RF to match gradient delay
        self._delay_rf()

        # add gradient ramps
        if add_ramps:
            self._add_ramps()

        # calculate modulation factor for slice offset
        self._calc_mod(use_girf)

    def verse_min_time(self, B_max, E_max=-1, add_ramps=True, use_girf=False):
        """
        Calculate minimum-time VERSE pulses from input rf pulse (Conolly 1988)
        Calculations are done on rf raster time (1us) and waveforms are interpolated in the end

        B_max: maximum RF amplitude [Hz]
        use_girf: Use girf for modulation function calculation
        """

        # energy constraint variables
        bmaxc = B_max
        bmaxh = B_max
        bmaxl = 0.0
        maxiter = 1000
        emaxratio=0.98
        niter = 0
        doneiter = False

        while not doneiter:
            niter += 1

            # calculate minimum time VERSE waveforms
            self._t_vs = self.rf_raster_time * np.maximum(abs(self.G)/self.max_grad, abs(self.rf)/bmaxc)
            self._t_vs[abs(self.rf)==0] = self.rf_raster_time
            self._alpha = np.divide(self.rf_raster_time, self._t_vs, out=np.ones_like(self._t_vs), where=self._t_vs>0)
            self.rf_vs = self.rf * self._alpha
            self.g_vs = self.G * self._alpha
            
            # adress possible slewrate overflow
            self._smooth_gradient()

            # interpolate to rf and gradient raster
            raster_ratio = round(self.grad_raster_time/self.rf_raster_time)
            vs_raster = np.cumsum(self._t_vs)
            rf_raster =  np.arange(0.5,self.N_vs+0.5) * self.rf_raster_time
            g_raster = np.arange(0.5,self.N_vs//raster_ratio+0.5) * self.grad_raster_time
            self.rf_vs = self._intp_wf(self.rf_vs, vs_raster, rf_raster)
            self.g_vs = self._intp_wf(self.g_vs, vs_raster, g_raster)
            
            if E_max > 0:
                benergy = self.calc_energy(use_input=False)
                if benergy > E_max:
                    bmaxh = bmaxc
                else:
                    bmaxl = bmaxc
                    if benergy/E_max > emaxratio:
                        doneiter = True
                    if niter==1:
                        print("Energy constraint fulfilled after 1st iteration. Exiting.\n")
                        bmaxl = 0
                        doneiter = True
                bmaxc = bmaxl + (bmaxh - bmaxl)/2.0

                if bmaxl >= bmaxh:
                    print("Warning:  Low-E limit > high-E limit.\n")
                    print(f"Exiting after {niter} iterations.\n")
                    doneiter = True
                if niter >= maxiter:
                    print("Warning:  Iteration limit reacheed.\n")
                    print(f"Exiting after {niter} iterations.\n")
                    doneiter = True
            else:
                doneiter = True

        print(f"VERSE pulse length: {1e3*self.N_vs*self.rf_raster_time} ms")

        # delay RF to match gradient delay
        self._delay_rf()

        # add gradient ramps
        if add_ramps:
            self._add_ramps()

        # calculate modulation factor for slice offset
        self._calc_mod(use_girf)

    def mintverse(self, B_max, E_max=-1, ramp_grad=False, use_girf=False):

        """
        Wraps mintverse C-code (Hargreaves, 2004)

        B_max: maximum RF amplitude [Hz]
        E_max: maximum RF energy (in units of dt*sum(B1^2), -1 to not constrain)
        ramp_grad: if True, ramp up/down gradient during RF pulse, otherwise ramps are added outside of RF pulse
                   Warning: should be False! ramping during RF doesnt seem to work properly
        use_girf: Use girf for modulation function calculation
        """

        import ctypes as ct
        import os
        mintv = ct.cdll.LoadLibrary(os.path.abspath("pulses/minverse/mintverse.so"))

        # zero-filling at start and end of rf/gradient to start/end them at 0
        rfr = self.rf.astype(np.complex64).real
        rfi = self.rf.astype(np.complex64).imag
        g = self.G
        rfr = np.concatenate(([0], rfr, [0]))
        rfi = np.concatenate(([0], rfi, [0]))
        if ramp_grad:
            g = np.concatenate(([0], g, [0]))
        else:
            g = np.concatenate((g[0], g, g[-1]))

        dt = np.ones_like(g) * self.rf_raster_time
        nout = np.int64([0])

        ob1r_data = (ct.c_double * 10000) ()
        ob1r = ct.pointer(ob1r_data)
        ob1r = ct.cast(ob1r, ct.POINTER(ct.POINTER(ct.c_double)))
        ob1i_data = (ct.c_double * 10000) ()
        ob1i = ct.pointer(ob1i_data)
        ob1i = ct.cast(ob1i, ct.POINTER(ct.POINTER(ct.c_double)))
        ogv_data = (ct.c_double * 10000) ()
        ogv = ct.pointer(ogv_data)
        ogv = ct.cast(ogv, ct.POINTER(ct.POINTER(ct.c_double)))

        raster_ratio = round(self.grad_raster_time/self.rf_raster_time)

        N_mod = 1
        while N_mod:
            res = mintv.mintverse(rfr.ctypes.data_as(ct.POINTER(ct.c_double)),
                    rfi.ctypes.data_as(ct.POINTER(ct.c_double)),
                    g.ctypes.data_as(ct.POINTER(ct.c_double)),
                    dt.ctypes.data_as(ct.POINTER(ct.c_double)),
                    len(dt),
                    ct.c_double(B_max),
                    ct.c_double(self.max_grad),
                    ct.c_double(self.max_slew),
                    ct.c_double(self.rf_raster_time),
                    ct.c_double(E_max),
                    nout.ctypes.data_as(ct.POINTER(ct.c_long)),
                    ob1r,
                    ob1i,
                    ogv)

            # ctypes to numpy
            self.N_vs = nout[0]
            _b1r = np.array([ob1r[0][_i] for _i in range(self.N_vs)])
            _b1i = np.array([ob1i[0][_i] for _i in range(self.N_vs)])
            self.rf_vs = _b1r + 1j*_b1i
            _g = np.array([ogv[0][_i] for _i in range(self.N_vs)])

            N_mod = self.N_vs % raster_ratio
            rfr = np.concatenate(([0], rfr, [0]))
            rfi = np.concatenate(([0], rfi, [0]))
            if ramp_grad:
                g = np.concatenate(([0,g[1]], g[1:-1], [g[-1],0]))
            else:
                g = np.concatenate((g[0], g, g[-1]))
            N_mod = 0


        print(f"VERSE pulse length: {1e3*self.N_vs*self.rf_raster_time} ms")

        # interpolate to gradient raster
        rf_raster =  np.arange(0.5,self.N_vs+0.5) * self.rf_raster_time
        g_raster = np.arange(0.5,self.N_vs//raster_ratio+0.5) * self.grad_raster_time
        self.g_vs = self._intp_wf(_g, rf_raster, g_raster)

        # delay RF to match gradient delay
        self._delay_rf()

        if not ramp_grad:
            self._add_ramps()

        # calculate modulation factor for slice offset
        self._calc_mod(use_girf)

    def mb_verse(self, band_sep, n_bands, phs_0_pt='None'):

        """ Modulate VERSE pulse for multiband imaging

        band_sep: band separation [m]
        n_bands: number of bands
        """
        from sigpy.mri.rf.multiband import mb_phs_tab

        if self.rf_vs is None:
            raise ValueError("Calculate VERSE pulse first.")

        if phs_0_pt != 'None':
            phs = mb_phs_tab(n_bands, phs_0_pt)
        else:
            phs = np.zeros(n_bands)

        # build multiband modulation function
        n = np.size(self.rf_vs)
        b = np.zeros(n, dtype='complex')
        for ii in range(0, n_bands):
            b += np.exp(1j* 2*np.pi * self.mod * band_sep * (ii - (n_bands - 1) / 2)) * np.exp(1j * phs[ii])

        self.rf_vs = self.rf_vs * b

    def _smooth_gradient(self):
        """Smooth verse gradient in case of a slewrate overflow (algorithm from Hargreaves' mintverse.c)
        """

        for k in range(len(self.g_vs)-1):
            # calculate averaged slewrate
            g_slew = 2 * (self.g_vs[k+1] - self.g_vs[k]) / (self._t_vs[k] + self._t_vs[k+1])
            if abs(g_slew) > self.max_slew:
                success = self._adjust_slew(k, self.max_slew)
                if not success:
                    break
        
        self._alpha = self.g_vs / self.G
        raster_ratio = round(self.grad_raster_time/self.rf_raster_time)
        self.N_vs = round(np.sum(self._t_vs)/self.grad_raster_time) * raster_ratio

    def _adjust_slew(self, ix, max_slew):
        # solve quadratic equation
        if self.g_vs[ix+1] >= self.g_vs[ix]:
            dth = self._t_vs[ix+1]
            dtl = self._t_vs[ix]
            gh = self.g_vs[ix+1]
            gl = self.g_vs[ix]
        else:
            dth = self._t_vs[ix]
            dtl = self._t_vs[ix+1]
            gh = self.g_vs[ix]
            gl = self.g_vs[ix+1]

        a = max_slew * dth
        b = max_slew * dtl + 2*gl
        c = -2 * gh

        d = b/2/a
        f = d*d - c/a

        if f < 0:
            print("Error: Quadratic roots are complex")
            stretch = -1
        else:
            g = np.sqrt(f)
            stretch = -d + g

        if self.g_vs[ix+1] >= self.g_vs[ix]:
            self.g_vs[ix+1] /= stretch
            self.rf_vs[ix+1] /= stretch
            self._t_vs[ix+1] *= stretch
        else:
            self.g_vs[ix] /= stretch
            self.rf_vs[ix] /= stretch
            self._t_vs[ix] *= stretch
            if stretch > 0 and ix > 0:
                g_slew = 2 * (self.g_vs[ix] - self.g_vs[ix-1]) / (self._t_vs[ix-1] + self._t_vs[ix])
                if abs(g_slew) > max_slew:
                    success = self._adjust_slew(ix-1, max_slew)
                    stretch = success + 0.5

        return stretch > 0
       
    def _intp_wf(self, wf, t_wf, t_to, mode='spline'):
        """
        interpolate waveforms by using cubic splines

        wf: waveform to interpolate
        t_wf: time raster of wf
        t_to: new time raster
        """

        if mode == 'spline':
            from scipy.interpolate import CubicSpline

            cs = CubicSpline(t_wf, wf)
            intp = cs(t_to)
        elif mode == 'linear':
            intp = np.interp(x=t_to, xp=t_wf, fp=wf)
        else:
            raise ValueError("Choose mode linear or spline.")

        return intp
    
    def _add_ramps(self):
        """
        add ramps to gradient
        """
        step = self.max_slew * self.grad_raster_time
        grad_edges = [self.g_vs[0], self.g_vs[-1]]
        ramps = [np.arange(0,grad_edges[0],step*np.sign(grad_edges[0]))]
        ramps.append(np.arange(0,grad_edges[1],step*np.sign(grad_edges[1]))[::-1])
        self.g_vs = np.concatenate((ramps[0],self.g_vs,ramps[1]))
        
        self.ramp_len = (len(ramps[0]), len(ramps[1])) # save ramp length for GIRF correction
      
    def _delay_rf(self):
        """
        Delay RF to match gradient delay
        """

        ndel = round(self.g_delay/self.rf_raster_time)
        self.rf_vs = np.concatenate([self.rf_vs,np.zeros(ndel)], axis=-1)
        self.rf_vs = np.roll(self.rf_vs, ndel, axis=-1)
        self.N_vs += ndel

    def _calc_mod(self, use_girf=False):
        """ Calculate modulation factor for slice offset
        """

        g_mod = self.g_vs
        if use_girf:
            g_mod = self._grad_pred(self.g_vs) # use GIRF? - gradient orientation unknown in advance
        if self.ramp_len[1] > 0:
            g_mod = g_mod[self.ramp_len[0]:-self.ramp_len[1]]

        rf_raster =  np.arange(0.5,len(self.rf_vs)+0.5) * self.rf_raster_time
        g_raster = np.arange(0.5,len(g_mod)+0.5) * self.grad_raster_time
        g_1us = self._intp_wf(g_mod, g_raster, rf_raster)
        self.mod = -1 * (np.cumsum(g_1us) - np.sum(g_1us)/2) * self.rf_raster_time # sign seems to be necessary for correct slice order

    def calc_energy(self, use_input=False):
        """ Calculates RF pulse energy

        use_input: Use input RF pulse instead of VERSE RF pulse
        """

        if use_input:
            rf = self.rf
        else:
            rf = self.rf_vs

        return np.sum(abs(rf**2)*self.rf_raster_time)

    ####
    # GIRF functions
    # Prediction is currently not really used as the orientation of the slice gradient is not known in advance
    # Instead, a simple gradient delay can be passed to the VERSE class
    # Also, the function _girf_corr is currently not correct
    ####

    def _grad_pred(self, grad_nom):
        """
        Predict gradient with GIRF - Does not work that well with gradients shaped similar as trapezoids
        Also, the orientation of the slice selection gradient is not known in advance, so GIRF correction is difficult

        grad_nom: nominal gradient on gradient raster time (10us)
        """

        # read girf
        girf = np.load("helper/girf_10us.npy")[-1,-1] # z-girf - only works for pure transverse orientation

        # predict gradient
        grad_smpl = len(grad_nom)
        grad_nom = np.concatenate((grad_nom, np.zeros(len(girf)-grad_smpl)))
        g_pred = fft(grad_nom)
        g_pred *= girf

        # cut out relevant part
        return ifft(g_pred)[:grad_smpl].real

    def _girf_corr(self):
        """ WIP/not working correctly
        This can not work as alpha is changing if the rf changes, see "verse_arbg.m" in Multiband-RF repo for correct implementation
        Also, the orientation of the slice selection gradient is not known in advance, so GIRF correction is difficult

        Correct RF pulse by predicting the actual gradient waveform with the GIRF
        """
        
        # predict gradient
        g_nom = self.g_vs.copy()
        g_pred = self._grad_pred(g_nom)

        # cut ramps
        g_nom = g_nom[self.ramp_len[0]:-self.ramp_len[1]]
        g_pred = g_pred[self.ramp_len[0]:-self.ramp_len[1]]

        # interpolate to rf raster
        rf_raster =  np.arange(0.5,len(self.rf_vs)+0.5) * self.rf_raster_time
        g_raster = np.arange(0.5,len(g_nom)+0.5) * self.grad_raster_time
        g_1us = self._intp_wf(g_nom, g_raster, rf_raster)
        g_pred_1us = self._intp_wf(g_pred, g_raster, rf_raster)
        
        # apply correction factor
        self.rf_vs *= g_pred_1us / g_1us

    def _grad_preemph(self):

        from helper.preEmphasis import preEmph
        
        girf = np.load("helper/girf_10us.npy")[-1,-1]
        grad_len = len(self.g_vs)
        grad = np.concatenate(([0], self.g_vs, [0]))
        grad_pre = preEmph(grad, girf, max_slewrate=self.max_slew)
        self.g_vs = grad_pre[1:grad_len+1]
