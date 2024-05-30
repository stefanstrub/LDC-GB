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

import lisaconstants as constants

from ldc.common.series import TimeSeries, FrequencySeries
from ldc.lisa.noise import simple_snr
import math
from ldc.lisa.orbits import AnalyticOrbits, Orbits
from ldc.waveform.waveform import GB_fdot
from ldc.waveform.fastGB import get_default_orbits, get_buffersize

import time

YEAR = constants.SIDEREALYEAR_J2000DAY*24*60*60
CLIGHT = constants.SPEED_OF_LIGHT

def computeXYZ(T, Gs, f0, fdot, fstar, ampl, q, tm, new=True):
    """ Compute TDI X, Y, Z from y_sr
    """
    f = f0 + fdot*tm
    omL = f/fstar
    SomL = np.sin(omL)
    fctr = np.exp(-1.j*omL)
    fctr2 = 4.0*omL*SomL*fctr/ampl

    ### I have factored out 1 - exp(1j*omL) and transformed to 
    ### fractional frequency: those are in fctr2
    ### I have rremoved Ampl to reduce dynamical range, will restore it later
    if new:
        Xsl =  Gs[:, 1, 0] - Gs[:, 2, 0] + (Gs[:, 0, 1] - Gs[:, 0, 2])*fctr
        Ysl =  Gs[:, 2, 1] - Gs[:, 0, 1] + (Gs[:, 1, 2] - Gs[:, 1, 0])*fctr
        Zsl =  Gs[:, 0, 2] - Gs[:, 1, 2] + (Gs[:, 2, 0] - Gs[:, 2, 1])*fctr
    else:
        Xsl =  Gs['21'] - Gs['31'] + (Gs['12'] - Gs['13'])*fctr
        Ysl =  Gs['32'] - Gs['12'] + (Gs['23'] - Gs['21'])*fctr
        Zsl =  Gs['13'] - Gs['23'] + (Gs['31'] - Gs['32'])*fctr
    XYZsl = fctr2*np.vstack([Xsl, Ysl, Zsl])
    XYZf_slow = ampl*np.fft.fft(XYZsl, axis=1)

    #for testing
    #Xtry = 4.0*(self.G21 - self.G31 + (self.G12 - self.G13)*fctr)/self.ampl

    M = XYZf_slow.shape[1] #len(XYZf_slow)
    XYZf = np.fft.fftshift(XYZf_slow, axes=1)
    f0 = (q - M/2) / T # freq = (q + np.arange(M) - M/2)/T
    return XYZf, f0

def construct_slow_part(T, arm_length, Ps, tm, f0, fdot, fstar, phi0, k, DP, DC, eplus, ecross, N=512):
    """
    Build linearly interpolated unit vectors of constellation arms.
    Matrix shape:
      * 3 satellites (vector start)
      * 3 satellites (vector end)
      * 3 coordinates
      * N points (linear interpolation of orbits)
    We only need to compute three vectors, since the matrix is skew-symmetric.
    TODO: Can we (or should we) even reduce it to a single vector ?

    """
    P = np.array(Ps)
    r = np.zeros(shape=(3, *P.shape))
    r[0, 1] = P[1] - P[0]
    r[0, 2] = P[2] - P[0]
    r[1, 2] = P[2] - P[1]
    r -= np.transpose(r, axes=(1, 0, 2, 3))
    r /= arm_length
    kdotr = np.dot(k, r)
    kdotP = np.transpose(np.dot(k, P), axes=(1, 0, 2)) / CLIGHT

    xi = tm - kdotP
    fi = f0 + fdot*xi
    # Ratio of true frequency to transfer frequency
    fonfs = fi / fstar 

    ### compute transfer f-n
    q = np.rint(f0 * T) # index of nearest Fourier bin
    df = 2.0 * np.pi * (q / T)
    om = 2.0 * np.pi * f0

    # Reduce to (3, 3, N).
    A = (eplus @ r * r * DP + ecross @ r * r * DC).sum(axis=2)

    ### Those are corrected values which match the time domain results.
    ## om*kdotP_i singed out for comparison with another code.
    argS =  phi0 + (om - df)*tm + np.pi*fdot*(xi**2)
    kdotP = om*kdotP - argS
    arg_ij = 0.5*fonfs * (1 + kdotr)
    G = 0.25*np.sinc(arg_ij/np.pi) * np.exp(-1.j*(arg_ij + kdotP)) * A
    return G, q


def construct_slow_part_old(T, arm_length, Ps, tm, f0, fdot, fstar, phi0, k, DP, DC, eplus, ecross, N=512):
    P1, P2, P3 = Ps
    r = dict()
    r['12'] = (P2 - P1)/arm_length ## [3xNt]
    r['13'] = (P3 - P1)/arm_length
    r['23'] = (P3 - P2)/arm_length
    r['31'] = -r['13']

    kdotr = dict()
    for ij in ['12', '13', '23']:
        kdotr[ij] = np.dot(k, r[ij]) ### should be size Nt
        kdotr[ij[-1]+ij[0]] = -kdotr[ij]

    kdotP = np.array([np.dot(k, P1),
                      np.dot(k, P2),
                      np.dot(k, P3)])
    kdotP /= CLIGHT
    
    Nt = len(tm)
    xi = tm - kdotP
    fi = f0 + fdot*xi
    fonfs = fi/fstar #Ratio of true frequency to transfer frequency

    ### compute transfer f-n
    q = np.rint(f0 * T) # index of nearest Fourier bin
    df = 2.*np.pi*(q/T)
    om = 2.0*np.pi*f0
    ### The expressions below are arg2_i with om*kR_i factored out

    A = dict()
    for ij in ['12', '23', '31']:
        aij = np.dot(eplus,r[ij])*r[ij]*DP+np.dot(ecross,r[ij])*r[ij]*DC
        A[ij] = aij.sum(axis=0)
    # These are wfm->TR + 1j*TI in c-code

    # arg2_1 = 2.0*np.pi*f0*xi[0] + phi0 - df*tm + np.pi*fdot*(xi[0]**2)
    # arg2_2 = 2.0*np.pi*f0*xi[1] + phi0 - df*tm + np.pi*fdot*(xi[1]**2)
    # arg2_3 = 2.0*np.pi*f0*xi[2] + phi0 - df*tm + np.pi*fdot*(xi[2]**2)

    ### These (y_sr) reproduce exactly the FastGB results 
    #self.y12 = 0.25*np.sin(arg12)/arg12 * np.exp(1.j*(arg12 + arg2_1)) * ( Dp12*self.DP + Dc12*self.DC )
    #self.y23 = 0.25*np.sin(arg23)/arg23 * np.exp(1.j*(arg23 + arg2_2)) * ( Dp23*self.DP + Dc23*self.DC )
    #self.y31 = 0.25*np.sin(arg31)/arg31 * np.exp(1.j*(arg31 + arg2_3)) * ( Dp31*self.DP + Dc31*self.DC )
    #self.y21 = 0.25*np.sin(arg21)/arg21 * np.exp(1.j*(arg21 + arg2_2)) * ( Dp12*self.DP + Dc12*self.DC )
    #self.y32 = 0.25*np.sin(arg32)/arg32 * np.exp(1.j*(arg32 + arg2_3)) * ( Dp23*self.DP + Dc23*self.DC )
    #self.y13 = 0.25*np.sin(arg13)/arg13 * np.exp(1.j*(arg13 + arg2_1)) * ( Dp31*self.DP + Dc31*self.DC )

    ### Those are corrected values which match the time domain results.
    ## om*kdotP_i singed out for comparison with another code.
    argS =  phi0 + (om - df)*tm + np.pi*fdot*(xi**2)
    kdotP = om*kdotP - argS
    Gs = dict()
    for ij, ij_sym, s in [('12', '12', 0), ('23', '23', 1), ('31', '31', 2),
                          ('21', '12', 1), ('32', '23', 2), ('13', '31', 0)]:
        arg_ij = 0.5*fonfs[s,:] * (1 + kdotr[ij])
        Gs[ij] = 0.25*np.sinc(arg_ij/np.pi) * np.exp(-1.j*(arg_ij + kdotP[s])) * A[ij_sym]
        
    ### Lines blow are extractions from another python code and from C-code
    # y = -0.5j*self.omL*A*sinc(args)*np.exp(-1.0j*(args + self.om*kq))
    # args = 0.5*self.omL*(1.0 - kn)
    # arg12 = 0.5*fonfs[0,:] * (1 + kdotr12)
    # arg2_1 = 2.0*np.pi*f0*xi[0] + phi0 - df*tm + np.pi*self.fdot*(xi[0]**2)  -> om*k.Ri
    # arg1 = 0.5*wfm->fonfs[i]*(1. + wfm->kdotr[i][j]);
    # arg2 =  PI*2*f0*wfm->xi[i] + phi0 - df*t;
    # sinc = 0.25*sin(arg1)/arg1;
    # tran1r = aevol*(wfm->dplus[i][j]*wfm->DPr + wfm->dcross[i][j]*wfm->DCr);
    # tran1i = aevol*(wfm->dplus[i][j]*wfm->DPi + wfm->dcross[i][j]*wfm->DCi);
    # tran2r = cos(arg1 + arg2);
    # tran2i = sin(arg1 + arg2);
    # wfm->TR[i][j] = sinc*(tran1r*tran2r - tran1i*tran2i);
    # wfm->TI[i][j] = sinc*(tran1r*tran2i + tran1i*tran2r);
    return Gs, q


cdef class pyGB:
    """ Galactic binary waveform fast generation. 
    """
    cdef public double arm_length
    cdef public double init_rotation
    cdef public double init_position
    cdef public long M,N
    cdef public double f0,fdot,ampl,theta,phi,psi,incl,phi0
    cdef public double T, delta_t
    cdef public int oversample
    cdef public int kmin
    cdef public object wfm
    cdef public object orbits
    
    def __cinit__(self, orbits=None, T=6.2914560e7, delta_t=15):
        """ Define C++ FastBinary dimensions and check that orbits are
        compatible.
        """
        if orbits is None or isinstance(orbits, dict):
            orbits = get_default_orbits()
        if not isinstance(orbits, AnalyticOrbits):
            raise TypeError('Fastbinary approximation requires analytic orbits') 
        self.orbits = orbits
        self.arm_length = self.orbits.arm_length
        self.T, self.delta_t = T, delta_t
        self.wfm = GB_fdot('fast', 'GB', 'GB_fdot')
        
    @property
    def citation(self):
        return '10.1103/physrevd.76.083006'
        
    def __reduce__(self):
        dorbits = dict({'orbit_type':'analytic',
                        'nominal_arm_length':self.orbits.arm_length, 
                        "initial_position": self.orbits.initial_position,
                        "initial_rotation": self.orbits.initial_rotation})
        return (self.__class__, (dorbits, self.T, self.delta_t))

    def info(self):
        return self.wfm.info()

    def buffersize(self, f0, ampl, oversample):
        """Get array dimension needed to compute TDI.
        """
        return get_buffersize(self.T, f0, ampl, oversample)

    def _parse_template(self, template):
        """Return source parameters from dictionary.
        """
        self.wfm.set_param(template)
        self.wfm.set_units()
        self.wfm.check_param()

        f0 = self.wfm.source_parameters["Frequency"]
        fdot = self.wfm.source_parameters["FrequencyDerivative"]
        theta = 0.5*np.pi-self.wfm.source_parameters['EclipticLatitude']
        phi = self.wfm.source_parameters['EclipticLongitude']
        ampl = self.wfm.source_parameters['Amplitude']
        incl = self.wfm.source_parameters['Inclination']
        psi = self.wfm.source_parameters['Polarization']
        phi0 = -self.wfm.source_parameters['InitialPhase']
        return [f0, fdot, ampl, theta, phi, psi, incl, phi0]

    def get_fd_tdixyz(self, radler=False, **kwargs):
        """ Return TDI X,Y,Z in freq. domain.
        
        f0 in Hz, fdot in Hz/s, ampl in strain,
        theta,phi,psi,incl,phi0 in rad.
        """
        if not radler:
            return self._get_fd_tdixyz_p(**kwargs)
        return self._get_fd_tdixyz_c(**kwargs)

    
    def _get_fd_tdixyz_p(self, template=None, f0=None, fdot=None, ampl=None,
                         theta=None, phi=None, psi=None, incl=None, phi0=None,
                         oversample=1, tdi2=False, freqs=None, xarray=True, new=True):
        """ Return TDI X,Y,Z in freq. domain.

        f0 in Hz, fdot in Hz/s, ampl in strain,
        theta,phi,psi,incl,phi0 in rad.
        """
        if template is not None:
            [f0, fdot, ampl, theta, phi, psi, incl, phi0] = self._parse_template(template)

        q0 = f0*self.T
        cosiota = np.cos(incl)
        fstar = CLIGHT/(self.arm_length*2*np.pi)
        cosps, sinps = np.cos(2.*psi), np.sin(2.*psi)
        Aplus = ampl*( 1.+ cosiota*cosiota)
        Across = -2.0*ampl*cosiota
        DP = Aplus*cosps - 1.0j*Across*sinps
        DC = -Aplus*sinps - 1.0j*Across*cosps

        sinth, costh = np.sin(theta), np.cos(theta)
        sinph, cosph = np.sin(phi), np.cos(phi)
        u = np.array([costh*cosph, costh*sinph, -sinth], ndmin=2)
        v = np.array([sinph, -cosph, 0], ndmin=2)
        k = np.array([-sinth*cosph, -sinth*sinph, -costh], ndmin=2)
        eplus = np.dot(v.T, v) - np.dot(u.T, u)
        ecross = np.dot(u.T, v) + np.dot(v.T, u)

        N = self.buffersize(f0, ampl, oversample) if freqs is None else len(freqs)
        tm = np.linspace(0, self.T, num=N, endpoint=False)
        Ps = [self.orbits.compute_position(1, tm),
              self.orbits.compute_position(2, tm),
              self.orbits.compute_position(3, tm)] # Pi = [3 x Nt] - coordinates vs time
        if new:
            Gs, q = construct_slow_part(self.T, self.arm_length, Ps, tm, f0, fdot, fstar,
                                        phi0, k, DP, DC, eplus, ecross, N=N)
        else:
            Gs, q = construct_slow_part_old(self.T, self.arm_length, Ps, tm, f0, fdot, fstar,
                                            phi0, k, DP, DC, eplus, ecross, N=N)
        Xf, f0 = computeXYZ(self.T, Gs, f0, fdot, fstar, ampl, q, tm, new=new)
        df = 1/self.T 
        kmin = int(np.round(f0/df))
        fctr = 0.5*self.T/N
        if tdi2:
            tdi2_factor = self.get_tdi2_factor(f0)
            fctr *= tdi2_factor
        Xf*= fctr
        if xarray:
            X, Y, Z = (FrequencySeries(Xf[0,:], df=df, kmin=kmin, t0=0, name="X"),
                       FrequencySeries(Xf[1,:], df=df, kmin=kmin, t0=0, name="Y"),
                       FrequencySeries(Xf[2,:], df=df, kmin=kmin, t0=0, name="Z"))
            return X, Y, Z
        return Xf, kmin

    def _get_fd_tdixyz_c(self, template=None, f0=None, fdot=None, ampl=None,
                         theta=None, phi=None, psi=None, incl=None, phi0=None,
                         oversample=1, simulator='synthlisa', tdi2=False,
                         freqs=None, xarray=True):
        """ Return TDI X,Y,Z in freq. domain.
        
        f0 in Hz, fdot in Hz/s, ampl in strain,
        theta,phi,psi,incl,phi0 in rad.
        """
        if template is not None:
            [f0, fdot, ampl, theta, phi, psi, incl, phi0] = self._parse_template(template)
        pars = [f0*self.T, np.cos(theta), phi, np.log(ampl),
                np.cos(incl), psi, phi0, fdot*self.T**2]

        if freqs is None :
            N = self.buffersize(f0,ampl,oversample)
        else : 
            N = len(freqs)
        M = N
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] xls = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] xsl = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] yls = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] ysl = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] zls = np.zeros(2*M)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] zsl = np.zeros(2*M)
        # TODO change to complex dtype

        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Cpars = np.array(pars)
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Opars = np.array([self.orbits.arm_length,
                                                                         self.orbits.initial_rotation,
                                                                         self.orbits.initial_position])
        Fast_GB_with_orbits(&Cpars[0], N, self.T, self.delta_t, &Opars[0],
                            &xls[0], &yls[0], &zls[0], &xsl[0], &ysl[0], &zsl[0],
                            len(pars))

        lout = [xsl, ysl, zsl] if simulator=="synthlisa" else [xls, yls, zls]
        fX,fY,fZ = [np.array(a[::2] + 1.j* a[1::2], dtype=np.complex128) for a in lout]
        kmin = int((f0*self.T) - M/2)
        df = 1.0/self.T
        #kmin = int(np.round(f0/df))
        if xarray:
            X, Y, Z = (FrequencySeries(fX/df, df=df, kmin=kmin, t0=0, name="X"),
                       FrequencySeries(fY/df, df=df, kmin=kmin, t0=0, name="Y"),
                       FrequencySeries(fZ/df, df=df, kmin=kmin, t0=0, name="Z"))
        else:
            X, Y, Z = fX/df, fY/df, fZ/df
        if tdi2:
            tdi2_factor = self.get_tdi2_factor(f0)
            X *= tdi2_factor
            Y *= tdi2_factor
            Z *= tdi2_factor
            
        return X,Y,Z

    def get_tdi2_factor(self, f0):
        omegaL = 2*np.pi*f0*(self.arm_length/CLIGHT)
        tdi2_factor = 2.j*np.sin(2*omegaL)*np.exp(-2j*omegaL)
        return tdi2_factor
    
    def get_td_tdixyz(self, **kwargs):
        """  Return TDI X,Y,Z in time domain.
        """
        fX, fY, fZ = self.get_fd_tdixyz(**kwargs)
        return (fX.ts.ifft(dt=self.delta_t),
                fY.ts.ifft(dt=self.delta_t),
                fZ.ts.ifft(dt=self.delta_t))
