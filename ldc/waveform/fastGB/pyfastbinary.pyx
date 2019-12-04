# distutils: sources = fastbinary.cc arrays.c
# distutils: language = c++
# cython: language_level=3

cdef extern from "fastbinary.hpp":
    cdef cppclass FastBinary:
        FastBinary(long, double, double) except +
        void response(double, double, double, double, double, double, double, double,
                      double*, double*, double*, double*, double*, double*, long, int)
        void XYZ(double, long, double*, double*, double*, double*, double*, double*)
        void setL(double)

        
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


year = constants.Nature.SIDEREALYEAR_J2000DAY*24*60*60

cdef class pyFastBinary:
    cdef FastBinary* FB
    cdef public double arm_length
    cdef public long M,N
    cdef public double f0,fdot,ampl,theta,phi,psi,incl,phi0
    cdef public double T, delta_t
    cdef public int oversample
    cdef public str method
    
    def __cinit__(self, orbits=None, template=None, T=6.2914560e7, delta_t=15, f0=None,
                  fdot=None, ampl=None, theta=None, phi=None,
                  psi=None, incl=None, phi0=None, oversample=1, method="mldc"):
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

        if template is not None:
            self._parse_template(template)
        else:
            self.f0, self.fdot, self.ampl = f0, fdot, ampl
            self.theta, self.phi = theta, phi
            self.psi, self.incl, self.phi0 = psi, incl, phi0
            
        M,N = self.buffersize(T, delta_t, method, oversample=oversample)
        self.FB = new FastBinary(long(N), T, delta_t)
        self.M, self.N = M, N
        self.T, self.delta_t = T, delta_t
        self.oversample = oversample
        self.method = method
        
    def __dealloc__(self):
        del self.FB

    def _set_mult(self, T, f):
        if T/year <= 1.0:
            mult = 1
        elif T/year <= 2.0:
            mult = 2
        elif T/year <= 4.0:
            mult = 4
        else:
            mult = 8

        # Doppler modulation, which is of order 1e-4 * f0
        # by contrast df = pi x 10^-8 Hz
        if f > 0.1:     # 631 bins/year
            N = 1024*mult
        elif f > 0.03:  # 189 bins/year
            N = 512*mult
        elif f > 0.01:  # 63 bins/year
            N = 256*mult
        elif f > 0.001: # 3 bins/year
            N = 64*mult
        else:
            N = 32*mult
        return mult, N
    
    def buffersize(self, T, dt, method, oversample):
        """
        """
        if method in ['legacy', 'mldc']:
            mult, N = self._set_mult(T, self.f0)

            # new logic for high-frequency chirping binaries (r1164 lisatools:Fast_Response.c)
            if method == 'mldc':
                try:
                    chirp = int(math.pow(2.0,math.ceil(math.log(self.fdot*T*T)/math.log(2.0))))
                except ValueError:
                    chirp = 0

                if chirp > N:
                    N = 2*chirp
                elif 4*chirp > N:
                    N = N * 2

            Acut = simple_snr(self.f0, self.ampl, years=T/year)

            # this is [2^(1+[log_2 SNR])];
            # remember that the response of the sinc is bound by 1/(delta bin #)
            # therefore we're using a number of bins
            # proportional to SNR, rounded to the next power of two
            # according to my tests, this is consistent with accepting +-0.4 in the detection SNR
            M = int(math.pow(2.0,1 + int(math.log(Acut)/math.log(2.0))));

            # corrected per Neil's 2009-07-28 changes to Fast_Response3.c
            # use the larger of the two buffer sizes, but not above 8192
            # -- further changed for new chirping logic (r1164)
            if method == 'mldc':
                M = N = max(M,N)
            else:
                M = N = min(8192,max(M,N))
        elif method=="lisasolve":  
            # LISA response bandwidth for Doppler v/c = 1e-4
            # on both sides, and frequency evolution
            deltaf = self.fdot * T + 2.0e-4 * self.f0
            # bins, rounded to next-highest power of two;
            # make sure we have enough for sidebands
            bins = 8 + deltaf*T
            N = int(2**math.ceil(math.log(bins,2)))

            # approximate source SNR
            f0 = self.f0 + 0.5 * self.fdot * T
            SNR = simple_snr(f0, self.ampl, years=T/year)

            # bandwidth for carrier frequency, off bin, worst case dx = 0.5 (accept 0.1 +- SNR)
            bins = max(1,2*2*SNR)
            M = int(2**math.ceil(math.log(bins,2)))
            M = N = max(M,N)
        else:
            self._check_method(method)
            
        M *= oversample; N *= oversample
        return M,N

    def _check_method(self, method):
        if not method in ["legacy", "mldc", "lisasolve"]:
            raise ValueError('Fastbinary method should be in [legacy, mldc, lisasolve]')

    def _parse_template(self, template):
        """ TODO: 
        - should be inherited from a general GB class
        - should check that keys exists
        - should also parse a vector ?
        
        """
        self.f0 = template["Frequency"]
        self.fdot = template["FrequencyDerivative"]
        self.theta = 0.5*np.pi-template['EclipticLatitude']
        self.phi = template['EclipticLongitude']
        self.ampl = template['Amplitude']
        self.incl = template['Inclination']
        self.psi = template['Polarization']
        self.phi0 = template['InitialPhase']
        
    def get_fd_tdixyz(self):
        """
        """

        cdef np.ndarray[np.double_t, ndim=1, mode="c"] xls = np.zeros(2*self.M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] xsl = np.zeros(2*self.M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] yls = np.zeros(2*self.M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] ysl = np.zeros(2*self.M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] zls = np.zeros(2*self.M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] zsl = np.zeros(2*self.M)
        # TODO change to complex dtype
        
        self.FB.response(self.f0, self.fdot, self.theta, self.phi, self.ampl, self.incl,
                         self.psi, self.phi0, 
                         &xls[0], &xsl[0], &yls[0], &ysl[0], &zls[0], &zsl[0], 2*self.M,
                         0 if self.method == 'legacy' else 1)

        fX,fY,fZ = [np.array(a[::2] + 1.j* a[1::2], dtype=np.complex128) for a in [xls, yls, zls]]
        # TODO convert to freq. array
        return fX,fY,fZ

    def get_td_tdixyz(self):
        fX, fY, fZ = self.get_fd_tdi()
        # todo compute inverse fft
        
