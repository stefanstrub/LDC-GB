# distutils: sources = GB.c LISA.c 
# distutils: language = c
# cython: language_level=3

cdef extern from "GB.h":
    void Fast_GB(double* , long, double, double,  double*, double*, double*, double*, double*, double*, int);
        

import numpy as np
cimport numpy as np
from ldc.common import constants
from ldc.lisa.noise import simple_snr
import math


# TODO:
# orbits: use ldc.orbits in fastbinary.cc
# arrays: use standard arrays
# parameters: check pycbc conventions, give info on expected units, getter/setter tools
# parameters for multiple sources and sum on TDI ?
# time domain
# noise
# output as frequency array


YEAR = constants.Nature.SIDEREALYEAR_J2000DAY*24*60*60

cdef class pyGB:
    cdef public double arm_length
    cdef public long M,N
    cdef public double f0,fdot,ampl,theta,phi,psi,incl,phi0
    cdef public double T, delta_t
    cdef public int oversample
    cdef public str method
    
    def __cinit__(self, orbits=None, T=6.2914560e7, delta_t=15):
        """ Define C++ FastBinary dimensions and check that orbits are compatible. 
        """
        if orbits is not None:
            if not isinstance(orbits, "AnalyticOrbits"):
                raise TypeError('Fastbinary approximation requires analytic orbits')
            else:
                self.arm_length = orbits.arm_length
                if orbits.initial_rotation !=0 or orbits.initial_position !=0:
                    raise ValueError('Fastbinary approximation requires null initial rotation and position')
        else:
            self.arm_length = 2.5e9
        self.T, self.delta_t = T, delta_t
        
    
    def buffersize(self, f0, ampl, method, oversample):
        """
        """
        Acut = simple_snr(f0,ampl,years=self.T/YEAR)
        mult = 8
        if((self.T/YEAR) <= 8.0): mult = 8
        if((self.T/YEAR) <= 4.0): mult = 4
        if((self.T/YEAR) <= 2.0): mult = 2
        if((self.T/YEAR) <= 1.0): mult = 1
        N = 32*mult
        if(f0 > 0.001): N = 64*mult
        if(f0 > 0.01):  N = 256*mult
        if(f0 > 0.03):  N = 512*mult
        if(f0 > 0.1):   N = 1024*mult

        M = int(math.pow(2.0,1 + int(np.log(Acut)/np.log(2.0))))
        if(M > 8192):
            M = 8192

        M = N = max(M,N)
        N *= oversample
        return(N)

    def _check_method(self, method):
        if not method in ["legacy", "mldc", "lisasolve"]:
            raise ValueError('Fastbinary method should be in [legacy, mldc, lisasolve]')

    def _parse_template(self, template):
        """ TODO: 
        - should be inherited from a general GB class
        - should check that keys exists
        - should also parse a vector ?
        
        """
        f0 = template["Frequency"]
        fdot = template["FrequencyDerivative"]
        theta = 0.5*np.pi-template['EclipticLatitude']
        phi = template['EclipticLongitude']
        ampl = template['Amplitude']
        incl = template['Inclination']
        psi = template['Polarization']
        phi0 = template['InitialPhase']
        return [f0, fdot, ampl, theta, phi, psi, incl, phi0]
        
    def get_fd_tdixyz(self, template=None, f0=None, fdot=None, ampl=None,
                      theta=None, phi=None, psi=None, incl=None, phi0=None,
                      method="mldc", oversample=1, simulator='synthlisa'):
        """
        """
        if template is not None:
            [f0, fdot, ampl, theta, phi, psi, incl, phi0] = self._parse_template(template)
        pars = [f0, fdot, ampl, theta, phi, psi, incl, phi0]

        N = self.buffersize(f0,ampl,method,oversample)
        M = N  
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] xls = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] xsl = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] yls = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] ysl = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] zls = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] zsl = np.zeros(2*M)
        # TODO change to complex dtype
        
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Cpars = np.array(pars)
        
        Fast_GB(&Cpars[0], N, self.T, self.delta_t,
                &xls[0], &yls[0], &zls[0], &xsl[0], &ysl[0], &zsl[0],
                len(pars))
            # self.f0, self.fdot, self.theta, self.phi, self.ampl, self.incl,
            #              self.psi, self.phi0, 
            #              &xls[0], &xsl[0], &yls[0], &ysl[0], &zls[0], &zsl[0], 2*self.M,
            #              0 if self.method == 'legacy' else 1)
        
        lout = [xsl, ysl, zsl] if simulator=="synthlisa" else [xls, yls, zls]
        fX,fY,fZ = [np.array(a[::2] + 1.j* a[1::2], dtype=np.complex128) for a in lout]
        kmin = int(f0*self.T) - M/2
        df = 1.0/self.T
        freq = np.linspace(kmin*df,(kmin + len(fX)-1)*df, len(fX))
        
        # TODO convert to freq. array
        return freq,fX,fY,fZ

    def get_td_tdixyz(self):
        fX, fY, fZ = self.get_fd_tdi()
        # todo compute inverse fft
        
