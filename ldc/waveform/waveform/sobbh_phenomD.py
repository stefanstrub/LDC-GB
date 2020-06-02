import numpy as np
import pyfftw
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from ldc.waveform.waveform.hphc import HpHc

class SOBBH_phenomD(HpHc):
    """ Compute waveforms h+ and hx of a stellar origin binary black hole. 

    Vectorized sources are not supported in this case. 
    """
    @property
    def dt(self):
        return self.source_parameters['Cadence']
    @property
    def redshift(self):
        return self.source_parameters['Redshift']
    @property
    def phi0(self):
        return self.source_parameters['InitialPhase']
    @property
    def incl(self):
        return self.source_parameters['Inclination']
    @property
    def m1s(self):
        return self.source_parameters['Mass1']
    @property
    def m2s(self):
        return self.source_parameters['Mass2']
    @property
    def chi1s(self):
        return self.source_parameters['Spin1']
    @property
    def chi2s(self):
        return self.source_parameters['Spin2']
    @property
    def Stheta1s(self):
        return self.source_parameters['PolarAngleOfSpin1']
    @property
    def Stheta2s(self):
        return self.source_parameters['PolarAngleOfSpin2']
    @property
    def DL(self):
        return self.source_parameters['Distance']*1e3 #Gpc -> Mpc
    @property
    def fstart(self):
        return self.source_parameters['InitialFrequency']
    @property
    def Tobs(self):
        return self.source_parameters['ObservationDuration']

    def precomputation(self):
        """ Load required parameters and convert them in expected units. """
        super().precomputation()
        self.a1 = np.cos(self.Stheta1s)*self.chi1s # For PhenomD we will use projections
        self.a2 = np.cos(self.Stheta2s)*self.chi2s
        


    def info(self):
        SOBBHunits = { "AzimuthalAngleOfSpin1": "Radian", 
                       "AzimuthalAngleOfSpin2": "Radian", 
                       "Distance": "Gpc", 
                       "EclipticLatitude": "Radian", 
                       "EclipticLongitude": "Radian", 
                       "Inclination": 'Radian', 
                       "InitialFrequency": 'Hz', 
                       "InitialPhase": 'Radian', 
                       "Mass1": 'SolarMass', 
                       "Mass2": 'SolarMass', 
                       "PolarAngleOfSpin1": "Radian", 
                       "PolarAngleOfSpin2": "Radian", 
                       "Polarization": "Radian", 
                       "Redshift": 'dimensionless', 
                       "Spin1": "MassSquared", 
                       "Spin2": "MassSquared",
                       'Cadence': 'Seconds',
                       'ObservationDuration': 'Seconds'}
        return SOBBHunits

    def check_param(self):
        for k in list(self.info().keys()):
            assert k in self.pnames
        assert self.units["EclipticLatitude"].lower() in ["radian", "rad", "r"]
        assert self.units["EclipticLongitude"].lower() in ["radian", "rad", "r"]
        assert self.units["PolarAngleOfSpin1"].lower() in ["radian", "rad", "r"]
        assert self.units["PolarAngleOfSpin2"].lower() in ["radian", "rad", "r"]
        assert self.units["Distance"].lower() in ["gpc"]
        assert self.units["Mass1"].lower() in ["solarmass"]
        assert self.units["Mass2"].lower() in ["solarmass"]



    def set_FD(self, fny=None, factor=2., maxf=1e-1, minf=1e-5):
        """ Set Fourier domain parameters
        """
        if fny is None:
            self.fny = 1.0/factor/self.dt      
        else:
            self.fny = fny

        self.maxf = maxf
        self.minf = minf


    def compute_hphc_td(self, t, source_parameters=None, approx_t=False):
        """ Return hp, hx for a time samples in t.

        Source parameters can be updated at the same time. 

        >>> GW = HpHc.type("my-sobbh", "SOBBH", "PhenomD")
        >>> hp,hc = GW.compute_hphc_td(np.arange(0,100,10), pSOBBH)
        >>> print(hp[0:3], hc[0:3] )
        [8.85674586e-22 1.48253757e-22 5.07077540e-22] [2.02727972e-21 2.16522345e-21 2.77731207e-21]
        """
        if source_parameters != None:
            self.set_param(source_parameters)

        # Check the approximant and call appropriate function
        if (self.approximant == 'PhenomD'):
            self.set_FD()
            tm, hpS, hcS = self.phenomD_SOBBH()
        else:
            raise NotImplementedError
        hpS, hcS = self.source2SSB(hpS, hcS) # Convert to SSB
        if approx_t:
            self.t, self.hp, self.hc = tm, hpS, hcS
        else:
            self._interp_hphc(tm, hpS, hcS, t, kind="spline") # interpolation
        return (self.hp, self.hc)

    def phenomD_SOBBH(self, phi_ref=0., tRef=0., tobs=0., nptmin=1000):
        """ Return hp,hc in the source frame. 

        """
        import pyFDresponse as FreqResp
        w_fr, w_amp, w_ph = FreqResp.GenerateResamplePhenomD(phi_ref, self.fstart,
                                                             self.m1s, self.m2s, self.a1, self.a2,
                                                             self.DL, self.incl,
                                                             minf=self.fstart, maxf=self.fny, 
                                                             settRefAtfRef=True, tRef=tRef,
                                                             tobs=tobs, nptmin=nptmin)

        # get rid of possibly huge constant in the phase before interpolating
        tfspline = spline(w_fr, 1/(2.*np.pi)*(w_ph-w_ph[0])).derivative() 
        tfvec = tfspline(w_fr)
        fspl = spline(tfvec, w_fr)
        fend = fspl(self.Tobs)

        df = 1.0/self.Tobs  # sampling rate to del_t
        Nf = int(self.fny/df+1)
        freq, hptilde, hctilde = FreqResp.ComputehphcFromAmpPhase([w_fr, w_amp, w_ph],
                                                                  np.arange(Nf)*df,
                                                                  self.incl, self.phi0)
        hptilde = np.conjugate(hptilde)
        hctilde = np.conjugate(hctilde)
        fbeg, fend = max(freq[0], w_fr[0]), freq[-1]# TODO fend=w_fr[0] ?? min(freq[-1], fend)
        ibeg, iend = int((fbeg - freq[0])/df), int((fend - freq[0])/df)
        hptilde[:ibeg],hctilde[:ibeg] = 0j, 0j
        hptilde[iend:],hctilde[iend:] = 0j, 0j
        del freq, w_fr, w_amp, w_ph

        #hp = np.fft.irfft(hptilde)*(1.0/self.dt)
        #hc = np.fft.irfft(hctilde)*(1.0/self.dt)
        hp = pyfftw.interfaces.numpy_fft.irfft(hptilde, n=2*(Nf-1))*(1.0/self.dt) # ifft
        hc = pyfftw.interfaces.numpy_fft.irfft(hctilde, n=2*(Nf-1))*(1.0/self.dt) # ifft
        tm = np.arange(len(hp))*self.dt
        return tm, hp, hc

    def generate_resample_phenomD(self, minf=1e-5, maxf=1., ):
        """ Generate PhenomD waveform and resample it in preparation for the FD response
        """
        # Sampling targets
        Delta_lnf = 0.02
        Delta_t = 1./24 # yrs
        Delta_f = 0.002 # Hz

        # Maximum frequency covered by PhenomD
        # NOTE: we take some small margin of error here. Before this was
        # interpreted as outside range by the C code and the last sample
        # was not computed (was given a phase of 0, which messed up the
        # interpolation afterwards)
        MfCUT_PhenomD = 0.2 - 1e-7

        # Check args
        if minf>=maxf:
            raise ValueError("Error in GenerateResamplePhenomDStartTime: incompatible minf and maxf.")
        if minf<=0 and tobs<=0:
            raise ValueError("Error in GenerateResamplePhenomDStartTime: both minf and tobs set to 0 and ignored, does not know where to start !")

        ## TODO : import from pyFDresponse

    


    
if __name__ == "__main__":
    import doctest
    import numpy as np

    pSOBBH = dict({"AzimuthalAngleOfSpin1": 0.0,# "Radian"), 
                   "AzimuthalAngleOfSpin2": 0.0,# "Radian"), 
                   "Distance": 0.8217407069275701,# "Gpc"), 
                   "EclipticLatitude": 0.23339632679489664,# "Radian"), 
                   "EclipticLongitude": 1.1798,# "Radian"), 
                   "Inclination": 1.1508,# 'Radian'), 
                   "InitialFrequency": 0.0074076, #'Hz'), 
                   "InitialPhase": 1.2622, #'Radian'), 
                   "Mass1": 31.033,# 'SolarMass'), 
                   "Mass2": 19.918,# 'SolarMass'), 
                   "PolarAngleOfSpin1": 2.7329, #"Radian"), 
                   "PolarAngleOfSpin2": 2.2947, #"Radian"), 
                   "Polarization": 3.7217, #"Radian"), 
                   "Redshift": 0.16454, #'unitless'), 
                   "Spin1": 0.4684,# "MassSquared"), 
                   "Spin2": 0.979,# "MassSquared"),
                   'Cadence': 5.,# 's')})
                   'ObservationDuration': 95.}) # 'Seconds'
    
    doctest.testmod()
