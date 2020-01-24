from abc import ABC, abstractmethod
from ldc.common import constants
import numpy as np
from scipy import interpolate
import os

C = constants.Nature
CLIGHT = C.VELOCITYOFLIGHT_CONSTANT_VACUUM
year = C.SIDEREALYEAR_J2000DAY*24*60*60


def simple_snr(f, h, i=None, years=1.0, noise_model='SciRDv1'):
    if not i:
        h0 = h * np.sqrt(16.0/5.0)    # rms average over inclinations
    else:
        h0 = h * np.sqrt((1 + np.cos(i)**2)**2 + (2.*np.cos(i))**2)
    sens = get_noise_model(noise_model,f).sensitivity(includewd=years)
    snr = h0 * np.sqrt(years * 365.25*24*3600) / np.sqrt(sens)
    return snr


def get_noise_model(model, frq=None, **kwargs):
    """Return the noise instance corresponding to model.
    
    model can be: "Proposal", "SciRDv1", "SciRDdeg1", "MRDv1", 
    "mldc", "newdrs", "LCESAcall"
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

class Noise():
    
    def __init__(self, frq, wd=0):
        self.freq = frq
        self.wd = wd

    def set_wdconfusion(self, duration):
        day = 24.*60*60
        if ((duration!=0 and duration < day/year) or (duration > 10.0)):
            raise NotImplementedError
        self.wd = duration
        
    @property
    def arm_length(self):
        return  2.5e9 # LISA's arm meters
 
    def WDconfusion(self, option="X"):
        """ option can be X or AE
        """
        duration = self.wd
        if duration==0:
            return 0
        lisaLT = self.arm_length/CLIGHT
        x = 2.0 * np.pi * lisaLT * self.freq
        t = 4.0 * x**2 * np.sin(x)**2
        Sg_sens = GalNoise(self.freq, duration*year).galshape()
        SgX = t*Sg_sens
        if option=="A" or option=="E":
            return 1.5*SgX
        elif option=="XY":
            return -0.5*SgX
        else:
            return SgX

    def to_file(self, filename):
        """ Save noise PSD to file
        """
        TDIn = np.rec.fromarrays([self.freq,
                                  self.psd(option="X"),
                                  self.psd(option="XY")], names=["freq", "X","XY"])
        np.save(filename, TDIn)

    def _XY2AET(self, X, XY, option='A'):
        if option in ["A", "E"]:
            return X-XY # = E
        elif option in ["T"]:
            return X+2*XY

    
    
        
class NumericNoise(Noise):
    """ PSD from file. 
    """
    def __init__(self, psdarray):
        self._psd = np.empty_like(psdarray)
        self.freq = psdarray["freq"]
        for k in ["X", "XY"]:
            self._psd[k] = psdarray[k]

    @staticmethod
    def from_file(filename):
        psd = np.load(filename)
        assert 'freq' in psd.dtype.names
        assert 'X' in psd.dtype.names
        assert 'XY' in psd.dtype.names
        N = NumericNoise(psd)
        return N
        
    def psd(self, freq=None, option='X'):
        if option in ['A', 'E', 'T']:
            XX = self.psd(freq, option="X", includewd=includewd)
            XY = self.psd(freq, option="XY", includewd=includewd)
            return self._XY2AET(XX, XY, option)

        if freq is not None:
            if np.isscalar(freq):
                freq = np.array([freq])
            S = np.interp(freq, self.freq, self._psd[option])
        if self.wd:
            S += self.WDconfusion(option=option)
        return S
    
class AnalyticNoise(Noise):
    """Analytic approximation of the two components of LISA noise:
    acceleration noise and optical metrology system (OMS) noise
    """
    
    def __init__(self, frq, model="SciRDv1", wd=0):

        Noise.__init__(self, frq, wd=wd)
        Sloc = (1.7e-12)**2 # m^2/Hz
        Ssci = (8.9e-12)**2 # m^2/Hz
        Soth = (2.e-12)**2  # m^2/Hz

        self.DSoms_d = {'Proposal':(10.e-12)**2,
                       'SciRDv1':(15.e-12)**2,
                       'MRDv1':(10.e-12)**2}  # m^2/Hz
        self.DSa_a = {'Proposal':(3.e-15)**2,
                     'SciRDv1':(3.e-15)**2,
                     'MRDv1':(2.4e-15)**2} # m^2/sec^4/Hz

        # Acceleration noise
        Sa_a = self.DSa_a[model] *(1.0 +(0.4e-3/frq)**2)*(1.0+(frq/8e-3)**4) # in acceleration
        self.Sa_d = Sa_a*(2.*np.pi*frq)**(-4.) # in displacement
        
        # In relative frequency unit
        Sa_nu = self.Sa_d*(2.0*np.pi*frq/CLIGHT)**2
        self.Spm =  Sa_nu
        
        # Optical Metrology System
        self.Soms_d = self.DSoms_d[model] * (1. + (2.e-3/frq)**4) # In displacement
        Soms_nu = self.Soms_d*(2.0*np.pi*frq/CLIGHT)**2 # In relative frequency unit
        self.Sop =  Soms_nu
        self.freq = frq
        self.model = model
        
    def relative_freq(self):
        return self.Spm, self.Sop

    def displacement(self):
        return self.Sa_d, self.Soms_d

    def sensitivity(self):#, includewd=0):
        ALL_m = np.sqrt( 4.0*self.Sa_d + self.Sop)
        AvResp = np.sqrt(5)         ## Average the antenna response
        Proj = 2./np.sqrt(3)     ## Projection effect

        ## Approximative transfert function
        f0 = 1.0/(2.0*(self.arm_length/CLIGHT))
        a = 0.41
        T  = np.sqrt(1+(self.freq/(a*f0))**2)
        Sens = (AvResp * Proj * T * ALL_m/self.arm_length)**2

        if (self.wd):
            day = 86400.0
            if ((self.wd < day/year) or (self.wd > 10.0)):
                raise NotImplementedError
            Sgal = GalNoise(self.freq, self.wd*year).galshape()
            Sens += Sgal
        return Sens

    def psd(self, freq=None, option="X"):#, includewd=0):
        """ option can be X, X2, XY, A, E, T
        """
        if option in ['A', 'E', 'T']:
            XX = self.psd(freq, option="X", includewd=includewd)
            XY = self.psd(freq, option="XY", includewd=includewd)
            return self._XY2AET(XX, XY, option)
            
        if self.wd and option in ["X2", "T"]:
            raise NotImplementedError

        if freq is not None:
            self.__init__(freq, model=self.model)
        
        lisaLT = self.arm_length/CLIGHT
        x = 2.0 * np.pi * lisaLT * self.freq
        if option=="X":
            S = 16.0 * np.sin(x)**2 * (2.0 * (1.0 + np.cos(x)**2) * self.Spm + self.Sop)
        elif option=="X2":
            S = 64.0 * np.sin(x)**2 * np.sin(2*x)**2 * self.Sop # TODO Check the acc. noise term
            S += 256.0 * (3 + np.cos(2*x)) * np.cos(x)**2 * np.sin(x)**4 * self.Spm
        elif option=="XY":
            S = -4.0 * np.sin(2*x) * np.sin(x) * (self.Sop + 4.0*self.Spm)
        elif option=="A" or option=="E":
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
            S += self.WDconfusion(option=option)
        return S

    
        
class SciRDdeg1Noise(AnalyticNoise):
    """
    """
    def __init__(self, frq, model="SciRDv1", wd=0):
        
        AnalyticNoise.__init__(self, frq, model, wd=wd)
        Sa_a = self.DSa_a['SciRDv1'] * (1.0+(0.4e-3/frq)**2) * \
               (1.0+(frq/8e-3)**4) * (1.0+(0.1e-3/frq)**2)
        self.Sa_d = Sa_a*(2.*np.pi*frq)**(-4.)
        Sa_nu = self.Sa_d*(2.0*np.pi*frq/CLIGHT)**2
        self.Spm =  Sa_nu


class SciRDdeg2Noise(AnalyticNoise):
    """
    """
    def __init__(self, frq, model="SciRDv1", wd=0):
        
        AnalyticNoise.__init__(self, frq, model, wd=wd)
        Sa_a = self.DSa_a['SciRDv1'] * (1.0+(0.4e-3/frq)**2) * \
               (1.0+(frq/8e-3)**4) * (1.0+(0.1e-3/frq)**4)
        self.Sa_d = Sa_a*(2.*np.pi*frq)**(-4.)
        Sa_nu = self.Sa_d*(2.0*np.pi*frq/CLIGHT)**2
        self.Spm =  Sa_nu

class MLDCNoise(AnalyticNoise):
    """
    """
    def __init__(self, f, wd=0):
        Noise.__init__(self, frq, wd=wd)
        self.Spm = 2.5e-48 * (1.0 + (f/1.0e-4)**-2) * f**(-2)
        defaultL = 16.6782
        lisaLT = self.arm_length/CLIGHT
        self.Sop = 1.8e-37 * (lisaLT/defaultL)**2 * f**2
        self.freq = f
        
    def displacement(self):
        raise NotImplementedError("")

class LCESAcallNoise(AnalyticNoise):
    """
    """
    def __init__(self, frq, wd=0):
        AnalyticNoise.__init__(self, frq, wd=wd)
        
        # Acceleration noise
        Sa_a = self.DSa_a['Proposal'] *(1.0 +(0.4e-3/frq)**2+(frq/9.3e-3)**4) # in acceleration
        self.Sa_d = Sa_a*(2.*np.pi*frq)**(-4.) # in displacement
        Sa_nu = self.Sa_d*(2.0*np.pi*frq/CLIGHT)**2 # in relative frequency unit
        self.Spm =  Sa_nu

        # Optical Metrology System
        self.Soms_d = self.DSoms_d['Proposal'] * (1. + (2.e-3/frq)**4) # in displacement
        Soms_nu = self.Soms_d*(2.0*np.pi*frq/CLIGHT)**2 # in relative frequency unit
        self.Sop =  Soms_nu
        self.freq = frq
        
class NewDRSNoise(AnalyticNoise):
    """
    """
    def __init__(self, f, wd=0):
        AnalyticNoise.__init__(self, frq, wd=wd)
        self.Spm = 6.00314e-48 * f**(-2) # 4.6e-15 m/s^2/sqrt(Hz)
        defaultL = 16.6782
        defaultD = 0.4
        defaultP = 1.0
       
        lisaD = 0.3  # TODO check it
        lisaP = 2.0  # TODO check it
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
    Here is an empirical fit for the Galactic confusion as a function of observation time
    We have used SNR>5.0 and SNR>7 as removal criteria.
    We use default to be 6 kinks, and used combined SNR. We also provide fits for 4 links
    @Tobs duration of observation (in seconds) should be between 3 months and 10 years
    """
    
    def __init__(self, freq, Tobs, SNR=7.0, links=6):

        self.freq = freq
        Tobs = Tobs/year

        if (SNR != 5.0 and SNR != 7.0):
            print ('We accept SNR to be 5.0 or 7.0', 'given', SNR)
            raise NotImplementedError

        self.Ampl, self.alpha, self.fr2, af1, bf1, afk, bfk = [ 1.28265531e-44, 1.62966700e+00,  4.81078093e-04, -2.23499956e-01,\
                                                                -2.70408439e+00, -3.60976122e-01, -2.37822436e+00]
        if (SNR == 5.0):
            self.Ampl, self.alpha, self.fr2, af1, bf1, afk, bfk =  [ 1.28265531e-44,  1.57926625e+00,  3.84466946e-04, -2.32198526e-01,\
                                                                    -2.77607250e+00, -2.78383832e-01, -2.51448127e+00]
        if (links == 4):
            self.Ampl, self.alpha, self.fr2, af1, bf1, afk, bfk =  [ 1.28265531e-44, 1.62966700e+00,  4.81078093e-04, -2.62516578e-01, \
                                                                    -2.60291132e+00, -4.20352978e-01, -2.23847341e+00]
            if (SNR == 5.0):
                self.Ampl, self.alpha, self.fr2, af1, bf1, afk, bfk =  [ 1.28265531e-44,  1.57926625e+00,  3.84466946e-04, -2.39762426e-01, -2.70133615e+00,\
                                                                        -3.49472046e-01, -2.37368814e+00]
        Tmin = 0.25
        Tmax = 10.0
        if (Tobs<Tmin or Tobs>Tmax):
            print ('Galaxy fit is valid between 3 months and 10 years, we do not extrapolate', Tobs, ' not in', Tmin, Tmax)
            raise NotImplementedError("")

        self.fr1 = 10.**(af1*np.log10(Tobs) + bf1)
        self.fknee = 10.**(afk*np.log10(Tobs) + bfk)


    def galshape(self):
        res = self.Ampl*np.exp(-(self.freq/self.fr1)**self.alpha) *\
              (self.freq**(-7./3.))*0.5*(1.0 + np.tanh(-(self.freq-self.fknee)/self.fr2))
        return res

        
        
