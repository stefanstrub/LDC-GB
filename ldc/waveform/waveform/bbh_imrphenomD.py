"""Compute waveforms h+ and hx for binary black holes using IMR
phenom D approximant."""

import numpy as np
import pyfftw
from scipy import signal

from ldc.common import constants
from ldc.common import tools
from ldc.waveform import imrphenomD
from ldc.waveform.waveform.hphc import HpHc

YRSID_SI = constants.Nature.SIDEREALYEAR_J2000DAY*24*60*60
MTsun = constants.Nature.SUN_GM/constants.Nature.VELOCITYOFLIGHT_CONSTANT_VACUUM**3

#pylint:disable=E1101
#pylint:disable=C0103

def butter_lowpass_filter(data, cutoff, fs, order=4, show=False):
    """ Return low pass filtered data 
    """
    #nyq = 0.5 * fs  # Nyquist Frequency
    #normal_cutoff = cutoff / nyq
    # Get the filter coefficients
    b, a = signal.butter(order, cutoff, btype='low', analog=False)
    y = signal.filtfilt(b, a, data)

    if show:
        import matplotlib.pyplot as plt
        w, h = signal.freqz(b, a)
        plt.figure()
        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    return y

class BBH_IMRPhenomD(HpHc):
    """ Compute waveforms h+ and hx of a massive black hole binary.

    Vectorized sources are not supported in this case.
    """
    parameter_map = {'redshift': 'Redshift',
                     'phi0': 'PhaseAtCoalescence',
                     'chi1s': 'Spin1',
                     'chi2s': 'Spin2',
                     'Stheta1s': 'PolarAngleOfSpin1',
                     'Stheta2s': 'PolarAngleOfSpin2',
                     'm1s': 'Mass1',
                     'm2s': 'Mass2',
                     'tc': 'CoalescenceTime',
                     'Tobs': 'ObservationDuration',
                     'dt': 'Cadence',
                     'DL': lambda p: p['Distance']*1e3} #Gpc -> Mpc}

    def precomputation(self):
        """ Load required parameters and convert them in expected units. """

        super().precomputation()
        p = self.source_parameters
        theL = p['InitialPolarAngleL']
        phiL = p['InitialAzimuthalAngleL']
        self.pol, self.incl = tools.aziPolAngleL2PsiIncl(self.eclLat, self.eclLon, theL, phiL)
        self.cos2Pol = np.cos(2.*self.pol)
        self.sin2Pol = np.sin(2.*self.pol)
        self.a1 = np.cos(self.Stheta1s)*self.chi1s # For PhenomD we will use projections
        self.a2 = np.cos(self.Stheta2s)*self.chi2s
        self.dist = self.DL*1e6*constants.Nature.PARSEC_METER
        self.set_FD()

    def info(self):
        """ Return default units
        """
        MBHBunits = {'EclipticLatitude':                 'Radian',\
                     'EclipticLongitude':                'Radian',\
                     'PolarAngleOfSpin1':                'Radian',\
                     'PolarAngleOfSpin2':                'Radian',\
                     'Spin1':                            'MassSquared',\
                     'Spin2':                            'MassSquared',\
                     'Mass1':                            'SolarMass',\
                     'Mass2':                            'SolarMass',\
                     'CoalescenceTime':                  'Second',\
                     'PhaseAtCoalescence':               'Radian',\
                     'InitialPolarAngleL':               'Radian',\
                     'InitialAzimuthalAngleL':           'Radian',\
                     'Cadence':                          'Seconds',\
                     'Redshift':                         'dimensionless',\
                     'Distance':                         'Gpc',
                     'ObservationDuration':              'Seconds'}
        return MBHBunits

    def check_param(self):
        """ Check parameters and their units
        """
        for k in list(self.info().keys()):
            assert k in self.pnames
        assert self.units["InitialPolarAngleL"].lower() in ["radian", "rad", "r"]
        assert self.units["EclipticLatitude"].lower() in ["radian", "rad", "r"]
        assert self.units["EclipticLongitude"].lower() in ["radian", "rad", "r"]
        assert self.units["InitialAzimuthalAngleL"].lower() in ["radian", "rad", "r"]
        assert self.units["PhaseAtCoalescence"].lower() in ["radian", "rad", "r"]
        assert self.units["PolarAngleOfSpin1"].lower() in ["radian", "rad", "r"]
        assert self.units["PolarAngleOfSpin2"].lower() in ["radian", "rad", "r"]
        assert self.units["Distance"].lower() in ["gpc"]
        assert self.units["Mass1"].lower() in ["solarmass"]
        assert self.units["Mass2"].lower() in ["solarmass"]
        assert self.units["CoalescenceTime"].lower() in ["s", "seconds", "second", "sec"]

        if not constants.check_cosmo(self.DL, self.redshift):
            print("DL and redshift do not match adopted cosmological model")
            raise ValueError


    def set_FD(self, df=None, factor=2., maxf_max=0.1):#0.8):
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

    def compute_hphc_td(self, t, source_parameters=None, approx_t=False, set_attr=True,
                        simbuffer=0., low_pass=False):
        """ Return hp, hx for a time samples in t.

        Source parameters can be updated at the same time.

        >>> GW = HpHc.type("my-mbhb", "MBHB", "IMRPhenomD")
        >>> hp,hc = GW.compute_hphc_td(np.arange(0,100,10), pMBHB)
        >>> print(hp[0:3], hc[0:3] )
        [ 1.01001832e-20  3.20629073e-20 -2.84964535e-21] [ 7.94954224e-21 -2.53257153e-20 -1.46996521e-20]
        """
        if source_parameters is not None:
            self.set_param(source_parameters)

        # Read param and convert in expected units

        # Check the approximant and call appropriate function
        if self.approximant == 'IMRPhenomD':
            tm, hpS, hcS = self.IMRPhenomD_MBHB(simbuffer=simbuffer)
        else:
            raise NotImplementedError

        if low_pass:
            cutoff = 1/(2*(t[1]-t[0])) #actual dt 5s -> cutoff at ~0.1Hz
            cutoff /= 1.25
            hpS = butter_lowpass_filter(hpS, cutoff, self.maxf)
            hcS = butter_lowpass_filter(hcS, cutoff, self.maxf)

        hpS, hcS = self.source2SSB(hpS, hcS) # Convert to SSB

        if approx_t:
            self.t, self.hp, self.hc = tm, hpS, hcS
        else:
            self._interp_hphc(tm, hpS, hcS, t, kind="spline") # interpolation
        return (self.hp, self.hc)

    def _IMRPhenomD_waveform(self, fRef=0., ph_ref=0., MfCUT_PhenomD=0.2-1e-7):
        """ Return ampl and phase in freq. domain
        """
        m1s_si = self.m1s*constants.Nature.SUN_MASS #  code needs masses in kg
        m2s_si = self.m2s*constants.Nature.SUN_MASS
        Mc = tools.mchirpofm1m2(self.m1s, self.m2s)
        f0 = tools.newtonianfoft(Mc, 2.0*self.Tobs/YRSID_SI)
        f0 = np.floor(f0/self.df)*self.df # redefine f0 to be an integer x df
        Ms = (self.m1s + self.m2s) * MTsun
        fmax = min(MfCUT_PhenomD/Ms, self.maxf)
        Nf = int((fmax-f0)/self.df)
        freq_all = np.arange(Nf)*self.df + f0 ### freq. array for waveform computation
        # Generate the signal
        wf_PhD_class = imrphenomD.pyIMRPhenomDh22AmpPhase(freq_all)
        frS, ampS, phS = wf_PhD_class.get_fd_waveform(ph_ref, fRef,
                                                      m1s_si, m2s_si,
                                                      self.a1, self.a2, self.dist)
        return frS, ampS, phS


    def IMRPhenomD_MBHB(self, simbuffer=0.):
        """ Return hp,hc in the source frame.

        Wrapper around MLDC IMRPhenomD model handling Fourier transforms.

        """
        frS, ampS, phS = self._IMRPhenomD_waveform()
        fmax = np.max(frS) # actual max Frequency
        if 0.5/fmax < self.dt:
            print("the recommended step size is smaller than user specified",
                  0.5/fmax, " < ", self.dt)
            raise ValueError

        factorp = 1/2.*(tools.spinWeightedSphericalHarmonic(-2, 2, 2, self.incl, self.phi0) + \
                          np.conj(tools.spinWeightedSphericalHarmonic(-2, 2, -2,
                                                                      self.incl, self.phi0)))
        factorc = 1j/2.*(tools.spinWeightedSphericalHarmonic(-2, 2, 2, self.incl, self.phi0) - \
                         np.conj(tools.spinWeightedSphericalHarmonic(-2, 2, -2,
                                                                     self.incl, self.phi0)))

        t0 = -self.tc - simbuffer # shift the data by tc
        phasetimeshift = -2.0*np.pi*frS*t0
        CommonPart = ampS * np.exp(1.0j*(phS+phasetimeshift)) ### 22 mode
        del ampS, phS, phasetimeshift

        # creating full arrays (from zero freq, to nyquist)
        fNy = 0.5/self.dt
        Nf = int(fNy//self.df)+1 #int(np.rint(fNy/self.df)+1)
        i_b = int(np.rint(frS[0]/self.df))
        i_e = i_b + len(frS)
        del frS
        hp_f = np.zeros(Nf, dtype='complex128')
        hc_f = np.zeros(Nf, dtype='complex128')
        hp_f[i_b:i_e] = np.conjugate(factorp * CommonPart)# to agree with fft conventions
        hc_f[i_b:i_e] = np.conjugate(factorc * CommonPart)
        Nfp = None #2*(Nf-1) if Nf%2!=0 else 2*Nf
        hpS = pyfftw.interfaces.numpy_fft.irfft(hp_f, n=Nfp)*(1.0/self.dt) # ifft
        hcS = pyfftw.interfaces.numpy_fft.irfft(hc_f, n=Nfp)*(1.0/self.dt) # ifft

        del hp_f, hc_f, CommonPart

        # taper and cut to Tobs
        Nt = len(hpS)
        tm = np.arange(Nt)*self.dt # corresponding time array
        i0 = np.argwhere((tm-self.Tobs) >= 0)[0][0]
        Ms = (self.m1s + self.m2s) * MTsun  # taper the end of the signal
        tap = self.tc + 500.0*Ms
        window = taper = 0.5*(1.0 - np.tanh(0.001*(tm  -  tap)))
        hpS, hcS = window*hpS, window*hcS
        return tm[:i0], hpS[:i0], hcS[:i0]

if __name__ == "__main__":
    import doctest
    pMBHB = dict({'EclipticLatitude': 0.312414, #"radian"),
                  'EclipticLongitude': -2.75291,# "radian"),
                  'CoalescenceTime': 28086000.0,# 's'),
                  'Distance':  9.14450149011798,# 'Gpc'),
                  'InitialAzimuthalAngleL': 3.9,# 'radian'),
                  'InitialPolarAngleL': 2.3535, #'radian'),
                  'Mass1': 132628.202,# "SolarMass"),
                  'Mass2': 30997.2481,# "SolarMass"),
                  'PhaseAtCoalescence':  3.8, #'Radian'),
                  'PolarAngleOfSpin1': 0.0,#'Radian'),
                  'PolarAngleOfSpin2': 0.0,#'Radian'),
                  'Redshift': 1.2687,# 'dimensionless'),
                  'Spin1': 0.9481998052314212, #'MassSquared'),
                  'Spin2': 0.9871324769575264, #'MassSquared'),
                  'Cadence': 5.,
                  'ObservationDuration':95.})#, 's') })
    doctest.testmod()
