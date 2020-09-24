# distutils: sources = GB.cc LISA.cc c_wrapper.cc
# distutils: language = c++
# cython: language_level=3

cdef extern from "GB.h":
    void Fast_GB(double* , long, double, double,
    	 	 double*, double*, double*, double*, double*, double*, int);
    void Fast_GB_with_orbits(double* , long, double, double, double*,
                             double*, double*, double*, double*, double*, double*, int);

        

import numpy as np
cimport numpy as np
from ldc.common import constants
from ldc.common.series import TimeSeries, FrequencySeries
from ldc.lisa.noise import simple_snr
import math
from ldc.lisa.orbits import AnalyticOrbits

# TODO:
# orbits: use ldc.orbits in lisa.c
# parameters: check pycbc conventions, give info on expected units, getter/setter tools



YEAR = constants.Nature.SIDEREALYEAR_J2000DAY*24*60*60

cdef class pyGB:
    cdef public double arm_length
    cdef public double init_rotation
    cdef public double init_position
    cdef public long M,N
    cdef public double f0,fdot,ampl,theta,phi,psi,incl,phi0
    cdef public double T, delta_t
    cdef public int oversample
    cdef public int kmin
    
    def __cinit__(self, orbits=None, T=6.2914560e7, delta_t=15):
        """ Define C++ FastBinary dimensions and check that orbits are
        compatible.
        """
        if orbits is not None:
            if not isinstance(orbits, AnalyticOrbits):
                raise TypeError('Fastbinary approximation requires analytic orbits')
            else:
                self.arm_length = orbits.arm_length
                self.init_rotation = orbits.initial_rotation
                self.init_position = orbits.initial_position
        else:
            self.arm_length = 2.5e9 # m
            self.init_rotation = 0 # rad
            self.init_position = 0 # rad
        self.T, self.delta_t = T, delta_t
    
    def buffersize(self, f0, ampl, oversample):
        """Get array dimension needed to compute TDI. 
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

    def _parse_template(self, template):
        """Return source parameters from dictionary.
        
        TODO: 
        - should be inherited from a general GB class
        - should check that keys exists
        - should also parse a vector ?
        - should get unit and check them
        """
        f0 = template["Frequency"]
        fdot = template["FrequencyDerivative"]
        theta = 0.5*np.pi-template['EclipticLatitude']
        phi = template['EclipticLongitude']
        ampl = template['Amplitude']
        incl = template['Inclination']
        psi = template['Polarization']
        phi0 = -template['InitialPhase']
        return [f0, fdot, ampl, theta, phi, psi, incl, phi0]
        
    def get_fd_tdixyz(self, template=None, f0=None, fdot=None, ampl=None,
                      theta=None, phi=None, psi=None, incl=None, phi0=None,
                      oversample=1, simulator='synthlisa'):
        """ Return TDI X,Y,Z in freq. domain. 
        
        f0 in Hz, fdot in Hz/s, ampl in strain, 
        theta,phi,psi,incl,phi0 in rad. 
        """
        if template is not None:
            [f0, fdot, ampl, theta, phi, psi, incl, phi0] = self._parse_template(template)
        pars = [f0*self.T, np.cos(theta), phi, np.log(ampl),
                np.cos(incl), psi, phi0, fdot*self.T**2]

        N = self.buffersize(f0,ampl,oversample)
        M = N  
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] xls = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] xsl = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] yls = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] ysl = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] zls = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] zsl = np.zeros(2*M)
        # TODO change to complex dtype
        
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Cpars = np.array(pars)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Opars = np.array([self.arm_length,
                                                                         self.init_rotation,
                                                                         self.init_position])
        Fast_GB_with_orbits(&Cpars[0], N, self.T, self.delta_t, &Opars[0],
                            &xls[0], &yls[0], &zls[0], &xsl[0], &ysl[0], &zsl[0],
                            len(pars))

        lout = [xsl, ysl, zsl] if simulator=="synthlisa" else [xls, yls, zls]
        fX,fY,fZ = [np.array(a[::2] + 1.j* a[1::2], dtype=np.complex128) for a in lout]
        kmin = int(int(f0*self.T) - M/2) 
        df = 1.0/self.T
        return (FrequencySeries(fX/df, df=df, kmin=kmin, t0=0, name="X"),
                FrequencySeries(fY/df, df=df, kmin=kmin, t0=0, name="Y"),
                FrequencySeries(fZ/df, df=df, kmin=kmin, t0=0, name="Z"))

    def get_td_tdixyz(self, **kwargs):
        """  Return TDI X,Y,Z in time domain. 
        """
        fX, fY, fZ = self.get_fd_tdixyz(**kwargs)
        return (fX.ts.ifft(dt=self.delta_t),
                fY.ts.ifft(dt=self.delta_t),
                fZ.ts.ifft(dt=self.delta_t))
    
