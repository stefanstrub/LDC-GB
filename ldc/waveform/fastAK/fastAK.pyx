# distutils: sources = EMRItemplate.cc, FreqAK_RAv2.cc
# distutils: language = c++
from libcpp cimport bool
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

from ldc.common.series import TimeSeries, FrequencySeries
from ldc.lisa.orbits import AnalyticOrbits

# cdef extern from "<complex.h>" namespace "std":
cdef extern from "EMRItemplate.h":
    cdef cppclass EMRItemplate:
          EMRItemplate() except +
          EMRItemplate(double, double, double, double, double) except +
          void SetPosition(double thetaS, double phiS, double thetaK, double phiK, double D)
          double M;
          double Mt;
          double m;
          double mt;
          double a;
          double lam;

          double e0;
          double nu0;
          double Phi0;
          double gamma0;
          double alpha0;
          double t0;  # instance of time when initial conditions are defined
          double fgam0;
          double falph0;

          double tPl;
          double e_pl;
          double nu_pl;
          double Phi_pl;
          double alpha_pl;
          double gamma_pl;

          double thS;   # co-latittude !!!
          double phS;
          double thK;
          double phK;
          double stS, stK, ctS, ctK; # cos and sin of thetaS and thetaK
          double cpS, spS, cpK, spK; # cos and sin of phiS, and phiK

          double dist; # distance in seconds
          double Ampl; # dimensionless amplitude: mu/D
          bool SkySet;

          double SNR;
          double LogL;

cdef extern from "FreqAK_RAv2.h":
    cdef cppclass FreqAK_RA:
          FreqAK_RA() except +
          FreqAK_RA(double, double, double, int, double*) except +
          void PhaseEv_RA(EMRItemplate& S, double t0, double dur, double tStart)
          void GetEvolution(double* t_phase, double* nus, double* eccs, double* phis,
                            double* alps, double* gams, double* fh);
          void ConstructWaveFreq_RA(EMRItemplate& S,
                                    double* Xf_r, double* Xf_im,
                                    double* Yf_r, double* Yf_im,
                                    double* Zf_r, double* Zf_im);
          void GetHarmonicsXYZ(EMRItemplate& S,
                               double** X_r, double** X_i,
                               double** Y_r, double** Y_i,
                               double** Z_r, double** Z_i,
                               double* phi, double* gam, double* alph,
                               double* tim, double** frqs, double** phase);
          void NewPhaseEv_RA(EMRItemplate& S, double t_ini, double dur, double tStart);
          void ComputeHAmpl(double ecc, int n, double* XaXb, double* Xa_Xb, double* Xc);
          void ConstructWaveFreq_RAv2(EMRItemplate& S,
                                      double* Xf_r, double* Xf_im,
                                      double* Yf_r, double* Yf_im,
                                      double* Zf_r, double* Zf_im);
          void FreeMemory();

cdef class pyFreqAK_RA:
    """ EMRI waveform fast generator. 
    """
    cdef FreqAK_RA fd_akra
    cdef EMRItemplate S
    cdef double maxDur
    cdef double timestep
    cdef double dtPhase
    cdef int Ext

    def __cinit__(self, delta_t=15, T=62914560, dtPhase=2048, extended=0, orbits=None):
        """

        delta_t: cadence for the time domain in s

        T: maximum duration of the observation, (stop waveform if
        plunge didn't occurs)
        
        dtPhase: cadence (internal) for evolving orbit and LISA
        motion, (usually 2048 sec).
        
        TODO: plug orbits and LISA constants
        """
        self.maxDur = T
        self.timestep = delta_t
        self.dtPhase = dtPhase
        self.Ext = extended

        if orbits is not None:
            if not isinstance(orbits, AnalyticOrbits):
                raise TypeError('Fastbinary approximation requires analytic orbits')
            else:
                arm_length = orbits.arm_length
                init_rotation = orbits.initial_rotation
                init_position = orbits.initial_position
        else:
            arm_length = 2.5e9 # m
            init_rotation = 0 # rad
            init_position = 0 # rad
        
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Opars = np.array([arm_length,
                                                                         init_rotation,
                                                                         init_position])
        self.fd_akra = FreqAK_RA(delta_t, T, dtPhase, extended, &Opars[0])

    @property
    def citation(self):
        return '10.1103/physrevd.69.082005'

    def _parse_template(self, pars):
        """
        pars is a dictionary
        TODO: set default units and check input units
    
        Distance in pc
        """
        if not "quadmom" in pars.keys() : pars["quadmom"] = 0.0
        self.S = EMRItemplate(pars['MBHMass'], pars['mu'], pars['spin'],
                             pars['lam'], pars['quadmom'])
        self.S.SetPosition(pars['thS'], pars['phS'], pars['thK'],
                           pars['phK'], pars['DL'])
        self.S.t0 = pars['t0']
        self.S.e0 = pars['e0']
        self.S.nu0 = pars['nu0']
        self.S.Phi0 = pars['phi0']
        self.S.alpha0 = pars['alph0']
        self.S.gamma0 = pars['gam0']

    def phase_ev_RA(self, t0, dur, t_start):
        """The signal will be generated tStart -> tStart+dur with IC given at
        t_start < t0 < t_start+dur
        
	Args:
            t0: moment of time for which intitial conditions are given
    	    dur: durarion of the signal from t_start
    	    t_start: how much back in time we should go from t0 -> t_start.
        """
        self.fd_akra.NewPhaseEv_RA(self.S, t0, dur, t_start)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_evolution(self):
        """ 
        """
        cdef int Nph = int(np.floor(self.maxDur / self.dtPhase ))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] t_phase = np.copy(np.zeros(Nph))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] nu = np.copy(np.zeros(Nph))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] ecc = np.copy(np.zeros(Nph))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] phi = np.copy(np.zeros(Nph))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] alp = np.copy(np.zeros(Nph))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] gam = np.copy(np.zeros(Nph))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] fH = np.copy(np.zeros(Nph))
        self.fd_akra.GetEvolution(&t_phase[0], &nu[0], &ecc[0], &phi[0],
                                  &alp[0], &gam[0], &fH[0])
        return (t_phase, nu, ecc, phi, alp, gam, fH)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_fd_tdixyz(self, template, t_start=0):
        """ Return TDI X,Y,Z in freq. domain.
        """

        self._parse_template(template)
        if t_start==0:
            self.phase_ev_RA(template["t0"], template["t0"], 0)
        
        cdef double fh_max = 0.5/self.timestep
        cdef int Npf = int(fh_max*self.maxDur)+1 # I allow to go beyond Nyquist freq. will be chopped off later

        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Xf_r = np.copy(np.zeros(Npf))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Xf_i = np.copy(np.zeros(Npf))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Yf_r = np.copy(np.zeros(Npf))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Yf_i = np.copy(np.zeros(Npf))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Zf_r = np.copy(np.zeros(Npf))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] Zf_i = np.copy(np.zeros(Npf))

        self.fd_akra.ConstructWaveFreq_RAv2(self.S, &Xf_r[0],  &Xf_i[0],
                                            &Yf_r[0], &Yf_i[0], &Zf_r[0], &Zf_i[0])
        Xf = Xf_r + 1.0j*Xf_i
        Yf = Yf_r + 1.0j*Yf_i
        Zf = Zf_r + 1.0j*Zf_i
        df = 1.0/self.maxDur
        Xf, Yf, Zf = (FrequencySeries(Xf, df=df, kmin=0, t0=0, name="X"),
                      FrequencySeries(Yf, df=df, kmin=0, t0=0, name="Y"),
                      FrequencySeries(Zf, df=df, kmin=0, t0=0, name="Z"))
        return (Xf, Yf, Zf)
        
    def get_td_tdixyz(self, **kwargs):
        """  Return TDI X,Y,Z in time domain.
        """
        fX, fY, fZ = self.get_fd_tdixyz(**kwargs)
        return (fX.ts.ifft(dt=self.timestep),
                fY.ts.ifft(dt=self.timestep),
                fZ.ts.ifft(dt=self.timestep))
    
    def destruct(self):
        self.fd_akra.FreeMemory()
