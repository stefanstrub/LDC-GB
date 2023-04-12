"""Compute waveforms h+ and hx for binary black holes using IMR
phenom D approximant."""

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline

import lisaconstants as constants

from ldc.common import tools
from ldc.common import constants as ldc_constants
from ldc.common.tools import butter_lowpass_filter, signal_inverse_spa
from ldc.common.series import FrequencySeries
from ldc.waveform import imrphenomD
from ldc.waveform.waveform.hphc import HpHc

#pylint:disable=E1101
#pylint:disable=C0103
#pylint:disable=W0201

YRSID_SI = constants.SIDEREALYEAR_J2000DAY*24*60*60
MTsun = constants.GM_SUN/constants.SPEED_OF_LIGHT**3

class BHB_IMRPhenomD(HpHc):
    """ Compute waveforms h+ and hx of a black hole binary either massive or of stellar origin.

    Vectorized sources are not supported in this case.
    """
    parameter_map = {'redshift':    'Redshift',
                     'phi0':        'PhaseAtCoalescence',
                     'm1s':         'Mass1',
                     'm2s':         'Mass2',
                     'chi1s':       'Spin1',
                     'chi2s':       'Spin2',
                     'Stheta1s':    'PolarAngleOfSpin1',
                     'Stheta2s':    'PolarAngleOfSpin2',
                     'tc':          'CoalescenceTime',
                     'Tobs':        'ObservationDuration',
                     'dt':          'Cadence',
                     'DL':          lambda p: p['Distance']} #Mpc -> Mpc}

    def precomputation(self):
        """ Load required parameters and convert them in expected units. """

        super().precomputation()
        p = self.source_parameters
        if self.source_type.lower() == "mbhb" :
            theL = p['InitialPolarAngleL']
            phiL = p['InitialAzimuthalAngleL']
            self.pol, self.incl = tools.aziPolAngleL2PsiIncl(self.eclLat,
                                                            self.eclLon,
                                                            theL,
                                                            phiL)
            self.cos2Pol = np.cos(2.*self.pol)
            self.sin2Pol = np.sin(2.*self.pol)

        if self.source_type.lower() == "sbbh":
            #Even though SOBHB are not defined by coalescence time,
            #this parameters allows to use compute_hphc_fd without any rewriting
            self.tc = 0.0
            self.incl = p['Inclination']
            self.phi0 = p['InitialPhase']
            self.f0 = p["InitialFrequency"]
            self.Stheta1s = 0.0
            self.Stheta2s = 0.0

        self.a1 = np.cos(self.Stheta1s)*self.chi1s # PhenomD projections
        self.a2 = np.cos(self.Stheta2s)*self.chi2s
        self.dist = self.DL*1e6*constants.PARSEC_METER
        self.set_FD()

    def info(self):
        """ Return default units
        """
        BHBunits= {
                    'Spin1':                    '1',\
                    'Spin2':                    '1',\
                    'Mass1':                    'Msun',\
                    'Mass2':                    'Msun',\
                    'EclipticLatitude':         'rad',\
                    'EclipticLongitude':        'rad',\
                    'Cadence':                  's',\
                    'ObservationDuration':      's',\
                    'Redshift':                 '1',\
                    'Distance':                 'Mpc',
                }

        if self.source_type.lower() == "mbhb" :
            BHBunits['CoalescenceTime'] = 's'
            BHBunits['PolarAngleOfSpin1'] = 'rad'
            BHBunits['PolarAngleOfSpin2'] = 'rad'
            BHBunits['PhaseAtCoalescence'] = 'rad'
            BHBunits['InitialPolarAngleL'] = 'rad'
            BHBunits['InitialAzimuthalAngleL'] = 'rad'

        if self.source_type.lower() == "sbbh":
            BHBunits["Inclination"] = 'rad'
            BHBunits["InitialPhase"] = 'rad'
            BHBunits["Polarization"] = 'rad'
            BHBunits["InitialFrequency"] = 'Hz'

        return BHBunits

    @property
    def citation(self):
        """ Return waveform citation info.
        """
        return imrphenomD._citation

    def check_param(self):
        """ Check parameters and their units
        """
        for k in list(self.info().keys()):
            # assert k in self.pnames
            try:
                assert k in self.pnames
            except AssertionError:
                raise AssertionError(f"Field {k} not found in input params")


        if self.source_type.lower() == "mbhb" :
            assert self.units["PhaseAtCoalescence"].lower() in ["radian", "rad", "r"]
            assert self.units["PolarAngleOfSpin1"].lower() in ["radian", "rad", "r"]
            assert self.units["PolarAngleOfSpin2"].lower() in ["radian", "rad", "r"]
            assert self.units["CoalescenceTime"].lower() in ["s", "seconds", "second", "sec"]
            assert self.units["InitialPolarAngleL"].lower() in ["radian", "rad", "r"]
            assert self.units["InitialAzimuthalAngleL"].lower() in ["radian", "rad", "r"]

        if self.source_type.lower() == "sbbh":
            assert self.units["Inclination"].lower() in ["radian", "rad", "r"]
            assert self.units["InitialPhase"].lower() in ["radian", "rad", "r"]
            assert self.units["Polarization"].lower() in ["radian", "rad", "r"]
            assert self.units["InitialFrequency"].lower() in ["hz"]

        assert self.units["EclipticLatitude"].lower() in ["radian", "rad", "r"]
        assert self.units["EclipticLongitude"].lower() in ["radian", "rad", "r"]
        assert self.units["Distance"].lower() in ["mpc"]
        assert self.units["Mass1"].lower() in ["solarmass", "msun", "solmass", "m_sun"]
        assert self.units["Mass2"].lower() in ["solarmass", "msun", "solmass", "m_sun"]
        assert self.units["Cadence"].lower() in ["s", "seconds", "second", "sec"]
        assert self.units["ObservationDuration"].lower() in ["s", "seconds", "second", "sec"]

        if not ldc_constants.check_cosmo(self.DL, self.redshift):
            print("DL and redshift do not match adopted cosmological model")
            raise ValueError

    def set_FD(self, df=None, factor=2., maxf_max=0.1):
        """ Set Fourier domain parameters
        """
        if df is None:
            # We use factor*Tobs to account for Gibbs oscillations and
            # will discard part of the signal.
            self.df = 1.0/(1.0*factor*self.Tobs)
        elif df is not None:
            self.df = df
        else:
            raise ValueError("Invalid parameter df or Tobs")
        self.maxf = 1.0/(1.0*self.dt) # define the max freq.
        if np.isscalar(self.maxf) and (self.maxf > maxf_max):
            self.maxf = maxf_max
        elif not np.isscalar(self.maxf):
            self.maxf[self.maxf > maxf_max] = maxf_max

    def MBHB_IMRPhenomD_waveform(self, fRef=0., ph_ref=0., MfCUT_PhenomD=0.2-1e-7):
        """ Return ampl and phase in freq. domain
        for a MBHB waveform
        """
        m1s_si = self.m1s*constants.SUN_MASS #  code needs masses in kg
        m2s_si = self.m2s*constants.SUN_MASS
        Mc = tools.mchirpofm1m2(self.m1s, self.m2s)
        f0 = tools.newtonianfoft(Mc, 2.0*self.Tobs/YRSID_SI)
        f0 = np.floor(f0/self.df)*self.df # redefine f0 to be an integer x df
        Ms = (self.m1s + self.m2s) * MTsun
        fmax = min(MfCUT_PhenomD/Ms, self.maxf)
        Nf = int((fmax-f0)/self.df)
        freq_all = np.arange(Nf)*self.df + f0 ### freq. array for waveform
        # Generate the signal
        PhD = imrphenomD.pyIMRPhenomDh22AmpPhase(freq_all)
        frS, ampS, phS = PhD.get_fd_waveform(ph_ref, fRef,m1s_si, m2s_si,
                                             self.a1, self.a2, self.dist)
        return frS, phS, ampS

    def SOBHB_frequency_grid(self, tf, frspl, ind, frS, nptmin) :
        """ Constructing the reduced freq. array for computing the response f-n.
        """
        # uniform in time
        Delta_t = 1.0/24.0*constants.SIDEREALYEAR_J2000DAY*86400.
        n_times = int(np.ceil(abs(tf[0] - tf[ind])/Delta_t))+1
        times_d_res = np.linspace(tf[0], tf[ind], n_times)

        freqs_d_res = frspl(times_d_res)
        freqs_d_res[0]= self.f0
        freqs_d_res[-1] = frS[-1]
        freqs_d_res =  freqs_d_res[freqs_d_res-frS[-1]<=0]

        # UNiform in freq.
        Delta_f = 5e-4
        npt_deltafsampling = int(np.ceil((frS[-1] - frS[0])/Delta_f))
        if npt_deltafsampling<=0:freq_rsdeltaf = np.array([])
        else: freq_rsdeltaf = frS[::int(len(frS)/npt_deltafsampling)]

        # uniform in log
        freq_rsln = np.logspace(np.log10(self.f0),np.log10(frS[-1]),num=nptmin)
        freq_tmp= np.concatenate((freqs_d_res, freq_rsdeltaf, freq_rsln[1:-1]))
        freq_tmp= np.unique(np.sort(freq_tmp))

        ind_fi= np.zeros(len(freq_tmp), dtype='int')
        freq_rs = np.zeros(len(freq_tmp))
        freq_rs[0]= freq_tmp[0]
        freq_rs[-1] = freq_tmp[-1]
        ind_fi[-1]= len(frS)-1
        for i in range(1, len(freq_tmp)-1):
            fi = freq_tmp[i]
            i_f = int(np.floor((fi-frS[0])*self.Tobs))
            freq_rs[i] = frS[i_f]
            ind_fi[i] = i_f

        freq_rs, ind_u = np.unique(freq_rs, return_index=True)
        ind_fi = ind_fi[ind_u]
        return freq_rs,ind_fi

    def SOBHB_IMRPhenomD_waveform(self, maxf=0.5, settRefAtfRef=True, tRef=0.0,
                                  nptmin=1000):
        """ Return ampl and phase in freq. domain for a SOBHB waveform
        """
        ## Generate the PhenomD phase and amplitude
        m1_SI = self.m1s*constants.SUN_MASS
        m2_SI = self.m2s*constants.SUN_MASS
        cnt = 2
        #### Generate PhenomD amplitude and phase on the fine and regular grid
        #### of points in freq for efficiency we guess the end frequency and
        #### extend it if not enough
        while cnt <= 50 :                   #TODO : Arbitrary, can be changed
            fmax= min(cnt*self.f0, maxf)
            if int(np.rint(fmax*self.Tobs)) == 1 :
                raise ValueError("The duration is not sufficient.")
            # freq_fromf0 = np.arange(self.f0,maxf,1./self.Tobs)
            freq_fromf0 = np.arange(self.f0,fmax,1./self.Tobs)
            # Generate the frequency, amplitude and phase of the waveform
            PhD = imrphenomD.pyIMRPhenomDh22AmpPhase(freq_fromf0)
            fr_a,amp_a,ph_a = PhD.get_fd_waveform(self.phi0, self.f0,
                                                    m1_SI, m2_SI,
                                                    self.a1, self.a2,
                                                    self.dist)

            # Spline the data
            tfspline= spline(fr_a, 1/(2.*np.pi)*ph_a).derivative()
            tf = tfspline(fr_a)
            tf -= tf[0]

            # Stopping the loop if max duration/frequency is reached
            if fmax == maxf:         #Check if max frequency is reached
                ind = np.argwhere(np.diff(tf)<0)
                ind = ind[0][0] if len(ind) else len(tf)-1
                frspl = spline(tf[:ind], fr_a[:ind])
                fr_end = fr_a[ind]
                break
            if tf[-1] > 2*self.Tobs: #Check if timerange is above duration
                ind = np.argwhere(tf > self.Tobs)[0][0]
                frspl = spline(tf[:ind], fr_a[:ind])
                fr_end = frspl(self.Tobs)
                break
            # TODO: Add default case so that fr_end is always defined somewhere.
            cnt += 1

        #Cut the end of data
        ampS,phS,frS= amp_a[fr_a<fr_end],ph_a[fr_a<fr_end],fr_a[fr_a<fr_end]

        #Get the frequency grid indexes
        freq_rs, ind_fi = self.SOBHB_frequency_grid(tf,frspl,ind,frS,nptmin)

        ### Shifting the to the reference
        if settRefAtfRef :
            fRef_inrange = min( max(frS[0], self.f0), frS[-1] )
            phS += 2.*np.pi*(frS - fRef_inrange)*(tRef - tfspline(fRef_inrange))
        return freq_rs,phS[ind_fi],ampS[ind_fi]

    def compute_hphc_fd(self, source_parameters=None, simbuffer=0., tapering=False):
        """ Return hp,hc in the source frame.

        Wrapper around MLDC IMRPhenomD model handling Fourier transforms.
        """
        #Ensure that all parameters are well written
        if source_parameters is not None:
            self.set_param(source_parameters)

        # Compute h22_fd
        if self.source_type.lower() == "mbhb" :
            frS, phS, ampS = self.MBHB_IMRPhenomD_waveform()
        if self.source_type.lower() == "sbbh" :
            frS, phS, ampS = self.SOBHB_IMRPhenomD_waveform()

        #Now let's convert to hp hc
        fmax = np.max(frS) # actual max Frequency

        if 0.5/fmax < self.dt:
            print("the recommended step size is smaller than user specified",
                  0.5/fmax, " < ", self.dt)
            raise ValueError

        # compute spherical harmonic coeffs associated to 22 mode
        spin_wsh_p = tools.spinWeightedSphericalHarmonic(-2, 2, 2,
                                                        self.incl, self.phi0)
        spin_wsh_m = tools.spinWeightedSphericalHarmonic(-2, 2, -2,
                                                        self.incl, self.phi0)

        K22p = 1/2.* ( spin_wsh_p + np.conj(spin_wsh_m) )
        K22c = 1j/2.*( spin_wsh_p - np.conj(spin_wsh_m) )

        #Phase and timeshift if necessary
        t0 = -self.tc - simbuffer # shift the data by tc
        phasetimeshift = -2.0*np.pi*frS*t0
        CommonPart = ampS * np.exp(1.0j * (phS))# + phasetimeshift)) ### 22 mode
        ExtraPart = np.exp(1.0j * (phasetimeshift))

        #Tapering before ifft
        if tapering:
            deltaminf, deltamaxf = 1e-5, 4e-3
            minf = tools.newtonianfoft(tools.mchirpofm1m2(self.m1s, self.m2s),
                                       self.Tobs/YRSID_SI)
            maxf = 5e-2
            w = tools.window_planck_vec(frS, minf, maxf, deltaminf, deltamaxf)
            CommonPart = w * CommonPart

        i_b = int(np.rint(frS[0]/self.df))
        del ampS, phS, phasetimeshift

        # creating full arrays (from zero freq, to nyquist)
        del frS
        hp_f = FrequencySeries(np.conjugate(K22p * CommonPart * ExtraPart),
                              df=self.df, kmin=i_b)
        hc_f = FrequencySeries(np.conjugate(K22c * CommonPart * ExtraPart),
                              df=self.df, kmin=i_b)

        return hp_f, hc_f

    def compute_hphc_td(self, t, source_parameters=None, approx_t=False,
                        set_attr=True, simbuffer=0., low_pass=False, tapering=True):
        """ Return hp, hx for a time samples in t.

        Source parameters can be updated at the same time.

        >>> MBHB = HpHc.type("my-mbhb", "MBHB", "IMRPhenomD")
        >>> hp_m,hc_m = MBHB.compute_hphc_td(np.arange(0,100,10), pMBHB)
        >>> SOBHB = HpHc.type("my-sobhb", "SBBH", "IMRPhenomD")
        >>> hp_s,hc_s = SOBHB.compute_hphc_td(np.arange(0,100,10),pSOBHB)
        """
        # Check the approximant and call appropriate function
        if 'IMRPhenom' in self.approximant and self.source_type.lower() == "mbhb" :
            hpS, hcS = self.compute_hphc_fd(source_parameters=source_parameters,
                                            simbuffer=simbuffer, tapering=tapering)

            # iFFT
            hpS = hpS.ts.ifft(dt=self.dt)
            hcS = hcS.ts.ifft(dt=self.dt)

        elif 'IMRPhenom' in self.approximant and self.source_type.lower() == "sbbh" :

            #Ensure that all parameters are well written
            if source_parameters is not None:
                self.set_param(source_parameters)

            # compute h22 in FD
            # fr, ph, amp = self.SOBHB_IMRPhenomD_waveform(nptmin=5000, tRef=tRef)
            fr, ph, amp = self.SOBHB_IMRPhenomD_waveform(nptmin=5000)

            # compute h22 in TD using SPA
            tspa, amp, ph = signal_inverse_spa(
                np.array([fr, amp, ph]).T, sign=1
            ).T

            # compute spherical harmonic coeffs associated to 22 mode
            spin_wsh_p = tools.spinWeightedSphericalHarmonic(-2, 2, 2,
                                                        self.incl, self.phi0)
            spin_wsh_m = tools.spinWeightedSphericalHarmonic(-2, 2, -2,
                                                        self.incl, self.phi0)

            K22p= 1/2. *( spin_wsh_p + np.conj(spin_wsh_m) )
            K22c= 1j/2.*( spin_wsh_p - np.conj(spin_wsh_m) )

            # interpolation function from inverse SPA timeline to wanted timeline
            spl_amp = spline(tspa, amp)
            spl_ph = spline(tspa, ph)

            ### NOK for phase
            ph_at0 = spl_ph(0.0)
            ph_at_t = spl_ph(t) - ph_at0

            # compute hp,hc for the wanted timeline
            hpS = spl_amp(t) * np.real( np.conjugate(K22p) * np.exp(-1.0j * ph_at_t) )
            hcS = spl_amp(t) * np.real( np.conjugate(K22c) * np.exp(-1.0j * ph_at_t) )

        else:
            raise NotImplementedError

        # From waveform source to SSB
        hpS, hcS = self.source2SSB(hpS, hcS) # Convert to SSB

        #Cut to Tobs
        tm = np.arange(len(hpS))*self.dt # corresponding time array
        # hcS, hpS, tm= hcS[tm<=self.Tobs], hpS[tm<=self.Tobs], tm[tm<=self.Tobs]

        if self.source_type.lower() == "mbhb" :
            # Taper
            #Caracteristic time for the merger
            Ms = (self.m1s + self.m2s) * MTsun
            #Getting rid of data long after the merger
            t_cutoff = min(tm[-1],self.tc + 1000*Ms)
            #Check that the data are long enough to perform the windowing
            if t_cutoff < tm[-1] :
                #Finding the maximum amplitude of both h+ and hx
                #and add 300 GM_sun/c**2*Ms to insure that enough
                # of the ringdown is taken into account
                tap = 0.5*(tm[tm>t_cutoff][np.argmax(np.array(hpS[tm>t_cutoff]))]
                           +tm[tm>t_cutoff][np.argmax(np.array(hpS[tm>t_cutoff]))])+300*Ms
                #Old window with new tap_time
                window = 0.5*(1.0 - np.tanh(0.001*(tm  -  tap)))
                hpS, hcS = window*hpS, window*hcS

        #Apply low pass filter
        if low_pass:
            cutoff = 1/(2*(t[1]-t[0])) #actual dt 5s -> cutoff at ~0.1Hz
            cutoff /= 1.25
            hpS = butter_lowpass_filter(hpS, cutoff, self.maxf)
            hcS = butter_lowpass_filter(hcS, cutoff, self.maxf)

        #Interpolation
        if approx_t:
            self.t, self.hp, self.hc = tm, hpS, hcS
            return self.hp, self.hc
        if self.source_type.lower() == "mbhb":
            self._interp_hphc(tm, hpS, hcS, t, kind="spline")
        elif self.source_type.lower() == "sbbh":
            self.t, self.hp, self.hc = t, hpS, hcS
        return self.hp, self.hc


if __name__ == "__main__":
    import doctest
    pMBHB = dict({'Mass1':                  132628.202,
                  'Mass2':                  30997.2481,
                  'Spin1':                  0.9481998052314212,
                  'Spin2':                  0.9871324769575264,
                  'EclipticLatitude':       0.312414,
                  'EclipticLongitude':      -2.75291,
                  'Redshift':               2.0,
                  'Distance':               15974.456786495544,
                  'Cadence':                5.,
                  'ObservationDuration':    95.,
                  'CoalescenceTime':        28086000.0,
                  'InitialAzimuthalAngleL': 3.9,
                  'InitialPolarAngleL':     2.3535,
                  'PhaseAtCoalescence':     3.8,
                  'PolarAngleOfSpin1':      0.0,
                  'PolarAngleOfSpin2':      0.0,
                  })
    pSOBHB = dict({"Mass1":                 50.,
                   "Spin1":                 0.0,
                   "Mass2":                 40.,
                   "Spin2":                 0.0,
                   "EclipticLatitude":      1.7,
                   "EclipticLongitude":     1.0471975511965976,
                   "Inclination":           1.0471975511965976,
                   "InitialFrequency":      1.2e-2,
                   "InitialPhase":          0.7,
                   "Polarization":          1.2,
                   "Redshift":              2.0,
                   "Distance":              15974.456786495544,
                   'Cadence':               5,
                   'ObservationDuration':   189348898.58127362/6})
    doctest.testmod()
