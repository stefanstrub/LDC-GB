import os
import numpy as np
from ldc.common import constants

C = constants.Nature
CLIGHT = C.VELOCITYOFLIGHT_CONSTANT_VACUUM
YEAR = C.SIDEREALYEAR_J2000DAY*24*60*60

known_noise_config = ["Proposal", "SciRDv1", "MRDv1", "mldc",
                      "SciRDdeg1", "SciRDdeg2", "newdrs"]


def simple_snr(freq, h, incl=None, years=1.0, noise_model='SciRDv1'):
    """Simple calculation using analytic noise approximation of LISA
    sensitivity.

    >>> simple_snr(1e-4, 1e-22)
    0.015476753564693699
    """
    if not incl:
        h0 = h * np.sqrt(16.0/5.0)    # rms average over inclinations
    else:
        h0 = h * np.sqrt((1 + np.cos(incl)**2)**2 + (2.*np.cos(incl))**2)
    sens = get_noise_model(noise_model, freq, wd=years).sensitivity()
    snr = h0 * np.sqrt(years * YEAR) / np.sqrt(sens)
    return snr


def get_noise_model(model, frq=None, **kwargs):
    """Return the noise instance corresponding to model.

    model can be: "Proposal", "SciRDv1", "SciRDdeg1", "MRDv1",
    "mldc", "newdrs", "LCESAcall"

    >>> N = get_noise_model("SciRDv1", np.logspace(-5, 0, 10000))
    """
    if model in ["Proposal", "SciRDv1", "MRDv1"]:
        NoiseClass = globals()["AnalyticNoise"]
        return NoiseClass(frq, model, **kwargs)
    elif model=="mldc":
        NoiseClass = globals()[model.upper()+"Noise"]
        return NoiseClass(frq, **kwargs)
    elif model=="SciRDdeg1":
        NoiseClass = globals()["SciRDdeg1Noise"]
        return NoiseClass(frq, **kwargs)
    elif model=="SciRDdeg2":
        NoiseClass = globals()["SciRDdeg2Noise"]
        return NoiseClass(frq, **kwargs)
    elif model=="newdrs":
        NoiseClass = globals()["NewRDSNoise"]
        return NoiseClass(frq, **kwargs)
    elif os.path.exists(model):
        NoiseClass = globals()["NumericNoise"]
        return NumericNoise.from_file(model)#NoiseClass(frq, model)
    else:
        NoiseClass = globals()[model+"Noise"]
        return NoiseClass(frq, **kwargs)

def _XY2AET(X, XY, option='A'):
    """ Convert TDI X and XY in TDI A, E, T.
    
    """
    if option in ["A", "E"]:
        return X-XY # = E
    elif option in ["T"]:
        return X+2*XY
    else:
        raise NotImplementedError

class Noise():

    def __init__(self, frq, wd=0):
        """Noise is defined for a given freq range.

        wd option add white dwarf confusion noise contribution, and is
        given in years.
        """
        self.freq = frq
        self.wd = wd

    def set_wdconfusion(self, duration):
        """ Set the white dwarf noise contribution in years.
        """
        day = 24.*60*60
        if ((duration!=0 and duration < day/YEAR) or (duration > 10.0)):
            raise NotImplementedError
        self.wd = duration

    @property
    def arm_length(self):
        return  2.5e9 # LISA's arm meters

    def wd_confusion(self, freq=None, option="X"):
        """ Compute the white dwarf confusion noise for PSD.

        option can be X, XY, A or E
        """
        if freq is None:
            freq = self.freq
        duration = self.wd
        if duration==0:
            return 0
        lisaLT = self.arm_length/CLIGHT
        x = 2.0 * np.pi * lisaLT * freq
        t = 4.0 * x**2 * np.sin(x)**2
        Sg_sens = GalNoise(freq, duration*YEAR).galshape()
        SgX = t*Sg_sens
        if option in ["A", "E"]:
            return 1.5*SgX
        elif option=="XY":
            return -0.5*SgX
        else:
            return SgX

    def to_file(self, filename):
        """Save noise PSD to file

        PSD are saved as numpy record array, and can be used later on
        with NumericNoise.

        >>> N = get_noise_model("SciRDv1", np.logspace(-5, 0, 10000))
        >>> fn = tmp_filename+".npy"
        >>> N.to_file(fn)
        >>> Nnum = get_noise_model(fn)
        >>> assert isinstance(Nnum, NumericNoise)
        """
        TDIn = np.rec.fromarrays([self.freq,
                                  self.psd(option="X"),
                                  self.psd(option="XY")], names=["freq", "X","XY"])
        np.save(filename, TDIn)


class NumericNoise(Noise):
    """ Noise PSD from file.
    """

    def __init__(self, psdarray): #pylint: disable=W0231
        """Set TDI noise and freq. attributes from numpy rec array.
        """
        self._psd = np.empty_like(psdarray)
        self.freq = psdarray["freq"]
        for k in ["X", "XY"]:
            self._psd[k] = psdarray[k]
        self.wd = 0

    @staticmethod
    def from_file(filename):
        """ Instantiate a Noise object from file.
        """
        psd = np.load(filename)
        assert 'freq' in psd.dtype.names
        assert 'X' in psd.dtype.names
        assert 'XY' in psd.dtype.names
        N = NumericNoise(psd)
        return N

    def psd(self, freq=None, option='X'):
        """Return noise PSD at given freq. or freq. range.

        Option can be X, XY, A, E, T.
        The input PSD read from file is interpolated at the given
        freq.
        """
        if option in ['A', 'E', 'T']:
            XX = self.psd(freq, option="X")
            XY = self.psd(freq, option="XY")
            return _XY2AET(XX, XY, option)

        if freq is not None:
            if np.isscalar(freq):
                freq = np.array([freq])
            S = np.interp(freq, self.freq, self._psd[option])
        if self.wd:
            S += self.wd_confusion(freq=freq, option=option)
        return S

class AnalyticNoise(Noise):
    """Analytic approximation of the two components of LISA noise:
    acceleration noise and optical metrology system (OMS) noise
    """

    def __init__(self, frq, model="SciRDv1", wd=0):
        """Set two components noise contributions wrt given model.
        """
        Noise.__init__(self, frq, wd=wd)

        self.DSoms_d = {'Proposal':(10.e-12)**2, 'SciRDv1':(15.e-12)**2,
                        'MRDv1':(10.e-12)**2}  # m^2/Hz
        self.DSa_a = {'Proposal':(3.e-15)**2, 'SciRDv1':(3.e-15)**2,
                      'MRDv1':(2.4e-15)**2} # m^2/sec^4/Hz

        # Acceleration noise
        Sa_a = self.DSa_a[model] * (1.0 +(0.4e-3/frq)**2) *\
               (1.0+(frq/8e-3)**4) # in acceleration
        self.Sa_d = Sa_a*(2.*np.pi*frq)**(-4.) # in displacement
        Sa_nu = self.Sa_d*(2.0*np.pi*frq/CLIGHT)**2 # in rel freq unit
        self.Spm =  Sa_nu

        # Optical Metrology System
        self.Soms_d = self.DSoms_d[model] * (1. + (2.e-3/frq)**4) # in displacement
        Soms_nu = self.Soms_d*(2.0*np.pi*frq/CLIGHT)**2 # in rel freq unit
        self.Sop =  Soms_nu

        self.freq = frq
        self.model = model

    def relative_freq(self):
        """ Return acceleration and OMS noise in relative freq unit
        """
        return self.Spm, self.Sop

    def displacement(self):
        """ Return acceleration and OMS noise in displacement
        """
        return self.Sa_d, self.Soms_d

    def sensitivity(self):
        """ Return sensitivity.

        >>> N = get_noise_model("SciRDv1", np.logspace(-5, 0, 5))
        >>> N.sensitivity()
        array([3.94460167e-27, 1.49307530e-34, 2.60680558e-40, 3.75087069e-41,
               9.95967435e-39])
        """
        all_m = np.sqrt(4.0*self.Sa_d + self.Sop)
        av_resp = np.sqrt(5) # average the antenna response
        proj = 2./np.sqrt(3) # projection effect

        ## Approximative transfert function
        f0 = 1.0 / (2.0* (self.arm_length/CLIGHT) )
        a = 0.41
        T = np.sqrt(1 + (self.freq/(a*f0))**2)
        sens = (av_resp * proj * T * all_m/self.arm_length)**2

        if self.wd:
            sens += GalNoise(self.freq, self.wd*YEAR).galshape()
        return sens

    def psd(self, freq=None, option="X"):#, includewd=0):
        """ Return noise PSD at given freq. or freq. range.

        Option can be X, X2, XY, A, E, T.

        >>> N = get_noise_model("SciRDv1", np.logspace(-5, 0, 5))
        >>> N.psd()
        array([7.13597299e-40, 2.76990908e-42, 9.52379492e-43, 1.92645601e-40,
               1.15359813e-36])
        """
        if option in ['A', 'E', 'T']:
            XX = self.psd(freq, option="X")
            XY = self.psd(freq, option="XY")
            return _XY2AET(XX, XY, option)

        if self.wd and option in ["X2", "T"]:
            raise NotImplementedError

        if freq is not None:
            self.__init__(freq, model=self.model, wd=self.wd)

        lisaLT = self.arm_length/CLIGHT
        x = 2.0 * np.pi * lisaLT * self.freq
        if option=="X":
            S = 16.0 * np.sin(x)**2 * (2.0 * (1.0 + np.cos(x)**2) * self.Spm + self.Sop)
        elif option=="X2":
            S = 64.0 * np.sin(x)**2 * np.sin(2*x)**2 * self.Sop
            S += 256.0 * (3 + np.cos(2*x)) * np.cos(x)**2 * np.sin(x)**4 * self.Spm
            # TODO Check the acc. noise term #pylint: disable=W0511
        elif option=="XY":
            S = -4.0 * np.sin(2*x) * np.sin(x) * (self.Sop + 4.0*self.Spm)
        elif option in ["A", "E"]:
            S = 8.0 * np.sin(x)**2 * (2.0 * self.Spm * (3.0 + 2.0*np.cos(x) +\
                                                        np.cos(2*x)) +\
                                      self.Sop * (2.0 + np.cos(x)))
        elif option=="T":
            S = 16.0 * self.Sop * (1.0 - np.cos(x)) * np.sin(x)**2 +\
                128.0 * self.Spm * np.sin(x)**2 * np.sin(0.5*x)**4
        else:
            print("PSD option should be in [X, X2, XY, A, E, T] (%s)"%option)
            return None
        if self.wd:
            S += self.wd_confusion(freq=freq, option=option)
        return S



class SciRDdeg1Noise(AnalyticNoise):
    """ First degration of SciRD 
    """
    def __init__(self, frq, model="SciRDv1", wd=0):

        AnalyticNoise.__init__(self, frq, model, wd=wd)
        Sa_a = self.DSa_a['SciRDv1'] * (1.0+(0.4e-3/frq)**2) * \
               (1.0+(frq/8e-3)**4) * (1.0+(0.1e-3/frq)**2)
        self.Sa_d = Sa_a*(2.*np.pi*frq)**(-4.)
        Sa_nu = self.Sa_d*(2.0*np.pi*frq/CLIGHT)**2
        self.Spm =  Sa_nu

class SciRDdeg2Noise(AnalyticNoise):
    """ Second degradation of SciRD
    """
    def __init__(self, frq, model="SciRDv1", wd=0):

        AnalyticNoise.__init__(self, frq, model, wd=wd)
        Sa_a = self.DSa_a['SciRDv1'] * (1.0+(0.4e-3/frq)**2) * \
               (1.0+(frq/8e-3)**4) * (1.0+(0.1e-3/frq)**4)
        self.Sa_d = Sa_a*(2.*np.pi*frq)**(-4.)
        Sa_nu = self.Sa_d*(2.0*np.pi*frq/CLIGHT)**2
        self.Spm =  Sa_nu

class MLDCNoise(AnalyticNoise):
    """ Noise model used in Radler
    """
    def __init__(self, f, wd=0):#pylint: disable=W0231
        Noise.__init__(self, f, wd=wd)#pylint: disable=W0233
        self.Spm = 2.5e-48 * (1.0 + (f/1.0e-4)**-2) * f**(-2)
        defaultL = 16.6782
        lisaLT = self.arm_length/CLIGHT
        self.Sop = 1.8e-37 * (lisaLT/defaultL)**2 * f**2
        self.freq = f

    def displacement(self):
        raise NotImplementedError("")

class LCESAcallNoise(AnalyticNoise):
    """ Noise model used for ESA call
    """
    def __init__(self, frq, wd=0):
        AnalyticNoise.__init__(self, frq, wd=wd)

        # Acceleration noise
        Sa_a = self.DSa_a['Proposal'] *\
               (1.0 +(0.4e-3/frq)**2+(frq/9.3e-3)**4) # in acceleration
        self.Sa_d = Sa_a*(2.*np.pi*frq)**(-4.) # in displacement
        Sa_nu = self.Sa_d*(2.0*np.pi*frq/CLIGHT)**2 # in relative frequency unit
        self.Spm =  Sa_nu

        # Optical Metrology System
        self.Soms_d = self.DSoms_d['Proposal'] * (1. + (2.e-3/frq)**4) # in displacement
        Soms_nu = self.Soms_d*(2.0*np.pi*frq/CLIGHT)**2 # in relative frequency unit
        self.Sop =  Soms_nu
        self.freq = frq

class NewDRSNoise(AnalyticNoise):
    """ ?
    """
    def __init__(self, f, wd=0):
        AnalyticNoise.__init__(self, f, wd=wd)
        self.Spm = 6.00314e-48 * f**(-2) # 4.6e-15 m/s^2/sqrt(Hz)
        defaultL = 16.6782
        defaultD = 0.4
        defaultP = 1.0

        lisaD = 0.3  
        lisaP = 2.0  
        lisaLT = self.arm_length/CLIGHT
        Sops = 6.15e-38 * (lisaLT/defaultL)**2 * (defaultD/lisaD)**4 *\
               (defaultP/lisaP) # 11.83 pm/sqrt(Hz)
        Sopo = 2.81e-38 # 8 pm/sqrt(Hz)
        self.Sop = (Sops + Sopo) * f**2
        self.freq = f

    def displacement(self):
        raise NotImplementedError("")


class GalNoise():
    """

    Provide an empirical fit for the galactic confusion as a function
    of observation time. We have used SNR>5.0 and SNR>7 as removal
    criteria.  We use default to be 6 links, and used combined SNR. We
    also provide fits for 4 links @Tobs duration of observation (in
    seconds) should be between 3 months and 10 years.

    """

    def __init__(self, freq, Tobs, snr=7.0, links=6):
        """ Set fitted values for galactic confusion.

        Tobs in given in seconds.
        """

        self.freq = freq
        Tobs = Tobs/YEAR

        if snr not in [5.0, 7.0]:
            print('We accept SNR to be 5 or 7', 'given', snr)
            raise NotImplementedError
        if links not in [6, 4]:
            print('We accept links to be 4 or 6', 'given', links)
            raise NotImplementedError

        L6_snr7 = [1.28265531e-44, 1.62966700e+00,  4.81078093e-04,
                   -2.23499956e-01, -2.70408439e+00, -3.60976122e-01, -2.37822436e+00]
        L6_snr5 = [1.28265531e-44,  1.57926625e+00,  3.84466946e-04,
                   -2.32198526e-01, -2.77607250e+00, -2.78383832e-01, -2.51448127e+00]
        L4_snr7 = [1.28265531e-44, 1.62966700e+00,  4.81078093e-04,
                   -2.62516578e-01, -2.60291132e+00, -4.20352978e-01, -2.23847341e+00]
        L4_snr5 = [1.28265531e-44,  1.57926625e+00,  3.84466946e-04,
                   -2.39762426e-01, -2.70133615e+00, -3.49472046e-01, -2.37368814e+00]
        if snr==5 and links==6:
            self.Ampl, self.alpha, self.fr2, af1, bf1, afk, bfk = L6_snr5
        elif snr==5 and links==4:
            self.Ampl, self.alpha, self.fr2, af1, bf1, afk, bfk = L4_snr5
        elif snr==7 and links==6:
            self.Ampl, self.alpha, self.fr2, af1, bf1, afk, bfk = L6_snr7
        elif snr==7 and links==4:
            self.Ampl, self.alpha, self.fr2, af1, bf1, afk, bfk = L4_snr7

        Tmin = 0.25
        Tmax = 10.0
        if (Tobs<Tmin or Tobs>Tmax):
            print('Galaxy fit is valid between 3 months and 10 years')
            print('we do not extrapolate', Tobs, ' not in', Tmin, Tmax)
            raise NotImplementedError("")

        self.fr1 = 10.**(af1*np.log10(Tobs) + bf1)
        self.fknee = 10.**(afk*np.log10(Tobs) + bfk)

    def galshape(self):
        """ Return galactic confusion contribution to noise PSD.
        """
        res = self.Ampl*np.exp(-(self.freq/self.fr1)**self.alpha) *\
              (self.freq**(-7./3.))*0.5*(1.0 + np.tanh(-(self.freq-self.fknee)/self.fr2))
        return res



if __name__ == "__main__":
    import doctest
    import tempfile
    tmp_filename = next(tempfile._get_candidate_names()) #pylint: disable=W0212
    doctest.testmod()
