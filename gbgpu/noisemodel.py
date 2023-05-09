# Noise curve taken from LDC used for Sangria and Spritz, "MRDv1"
from abc import ABC, abstractmethod
import numpy as cp
CLIGHT = 299792458.
arm_length = 2.5e9

class AnalyticNoise(ABC):
    """Analytic approximation of the two components of LISA noise:
    acceleration noise and optical metrology system (OMS) noise
    """

    def __init__(self, frq):
        """Set two components noise contributions wrt given model.
        """
        super().__init__()

        self.DSoms_d = (10.e-12)**2  # m^2/Hz
        self.DSa_a = (2.4e-15)**2 # m^2/sec^4/Hz
        self.freq = frq

        # Acceleration noise
        Sa_a = self.DSa_a * (1.0 +(0.4e-3/frq)**2) * (1.0+(frq/8e-3)**4) # in acceleration
        self.Sa_d = Sa_a*(2.*cp.pi*frq)**(-4.) # in displacement
        Sa_nu = self.Sa_d*(2.0*cp.pi*frq/CLIGHT)**2 # in rel freq unit
        self.Spm =  Sa_nu

        # Optical Metrology System
        self.Soms_d = self.DSoms_d * (1. + (2.e-3/frq)**4) # in displacement
        Soms_nu = self.Soms_d*(2.0*cp.pi*frq/CLIGHT)**2 # in rel freq unit
        self.Sop =  Soms_nu

    def psd(self, option="A", tdi2=False):#, includewd=0):
        """ Return noise PSD at given freq. or freq. range.

        Option can be X, A, E, T.

        """
        lisaLT = arm_length/CLIGHT
        x = 2.0 * cp.pi * lisaLT * self.freq #cp.array(self.freq)

        if option=="X":
            S = 16.0 * cp.sin(x)**2 * (2.0 * (1.0 + cp.cos(x)**2) * self.Spm + self.Sop)
        elif option in ["A", "E"]:
            S = 8.0 * cp.sin(x)**2  * (2.0 * self.Spm * (3.0 + 2.0*cp.cos(x) + cp.cos(2.0*x)) + self.Sop * (2.0 + cp.cos(x)))
        elif option=="T":
            S = 16.0 * self.Sop * (1.0 - cp.cos(x)) * cp.sin(x)**2 + 128.0 * self.Spm * cp.sin(x)**2 * cp.sin(0.5*x)**4
        else:
            print("PSD option should be in [X, A, E, T] (%s)"%option)
            return None
        if tdi2:
            factor_tdi2 = 4 * cp.sin(2 * x)**2
            S *= factor_tdi2
        return S

