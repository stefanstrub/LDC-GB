"""Compute waveforms h+ and hx for extreme mass ratio inspiral (EMRI) systems using
FastSchwarschildEccentricFlux (FSEF) model.
"""
import numpy as np
from astropy import units as un

from lisaconstants import SIDEREALYEAR_J2000DAY
from ldc.common import tools
from ldc.waveform.waveform.hphc import HpHc

#pylint:disable=E1101
#pylint:disable=C0103
#pylint:disable=C0303
#pylint:disable=W0201


class EMRI_FastSchwarschildEccentricFlux(HpHc):
    """ Compute waveforms h+ and hx of an extreme mass ratio inspiral (EMRI) system

    Vectorized sources are not supported in this case.
    """
    
    def __init__(self, *args):
        """Import few package and init a FastSchwarzschildEccentricFlux
        container.

        """
        from few.waveform import FastSchwarzschildEccentricFlux
        super().__init__(*args)
        self.few = FastSchwarzschildEccentricFlux()

    def precomputation(self):
        """ Load required parameters and convert them in expected units. 
        
        .. note :: PhiPhi0 = - phi_0 (InitialPhase)
        """
        super().precomputation()
        p = self.source_parameters
        
        self.theL = p['InitialPolarAngleL']
        self.phiL = p['InitialAzimuthalAngleL']
        self.pol, self.incl = tools.aziPolAngleL2PsiIncl(self.eclLat, self.eclLon,
                                                         self.theL, self.phiL)
        self.cos2Pol = np.cos(2.*self.pol)
        self.sin2Pol = np.sin(2.*self.pol)
        
        self.M = p['MassMBHB']
        self.mu = p["MassSOBHB"]
        self.p0 = p["InitialSemiLatusRect"]
        self.e0 = p["InitialEccentricity"]
        self.phi_r0 = p['InitialPhasePhiR']
        self.phi_phi0 = p['InitialPhase']
        self.dist = p["Distance"]
        self.dt = p["Cadence"]
        self.T = p["ObservationDuration"]
        self.eps = p["eps"]

    @property
    def citation(self):
        return 'arxiv.org/2104.04582;arxiv.org/2008.06071'
        
    def info(self):
        """ Return default units
        """
        EMRIfsefunits = {
            "MassMBHB":                             'Msun', # M
            "MassSOBHB":                            'Msun', # mu
            "InitialSemiLatusRect":                  '1',    # p0
            "InitialEccentricity":                  '1',    # e0
            'InitialPolarAngleL':                   'rad',
            'InitialAzimuthalAngleL':               'rad',
            'InitialPhasePhiR':                     'rad',
            'InitialPhase':                         'rad',
            "Distance":                             'Mpc',  # dist
            "Cadence":                              's',    # dt
            "ObservationDuration":                  's', # T
            "eps":                                  '1',    # eps
            "EclipticLatitude":                     'rad',
            "EclipticLongitude":                    'rad',
        }
        return EMRIfsefunits
    
    def check_param(self):
        """ Check parameters and their units.
        
        ..note :: parameter values are checked within the few lib.
        """
        for k in list(self.info().keys()):
            # assert k in self.pnames
            try:
                assert k in self.pnames
            except AssertionError:
                print(f"{k} param not present")
        assert self.units["MassMBHB"].lower() in ["msun"]
        assert self.units["MassSOBHB"].lower() in ["msun"]
        assert self.units["EclipticLatitude"].lower() in ["radian", "rad", "r"]
        assert self.units["EclipticLongitude"].lower() in ["radian", "rad", "r"]
        assert self.units["InitialPhasePhiR"].lower() in ["radian", "rad", "r"]
        assert self.units["InitialPhase"].lower() in ["radian", "rad", "r"]
        assert self.units["InitialPolarAngleL"].lower() in ["radian", "rad", "r"]
        assert self.units["InitialAzimuthalAngleL"].lower() in ["radian", "rad", "r"]
        assert self.units["Distance"].lower() in ["mpc"]
        assert self.units["ObservationDuration"].lower() in ["s"]

    def compute_hphc_td(self, t, source_parameters=None, approx_t=False, set_attr=False):
        """ Return hp, hx for a time samples in t.
        Source parameters can be updated at the same time.
        
        warning :: inconsistent notation with SOBHB: hpS for source frame and hp for obs frame

        :param t: time vector array
        :type t: numpy array or list
        :param source_parameters: source physical parameters
        :type source_parameters: dict
        :param approx_t: whether or not to use spline interpolation
        :type approx_t: boolean
        :return: hplus and hcross GW polarisations
        :rtype: (numpy.ndarray, numpy.ndarray) 
        
        """
        if source_parameters is not None:
            self.set_param(source_parameters)
            
        # convert from Mpc to Gpc (needed by FEW)
        dist = self.dist * 1e-3
        # convert from s to yrs (needed by FEW)
        T = self.T / (SIDEREALYEAR_J2000DAY*24*3600)

        # compute hplus and hcross in 
        # the source frame.
        wave = self.few(
            self.M, self.mu, self.p0, self.e0, 
            self.theL, self.phiL, 
            dist=dist,
            Phi_r0=self.phi_r0,
            Phi_phi0=self.phi_phi0,
            T=T, dt=self.dt
        )
        # wave is a complex vector h(t) = hp(t)-i*hc(t)
        # cf. FEW doc
        hpS, hcS = wave.real, - wave.imag

        # convert to SSB
        hp, hc = self.source2SSB(hpS, hcS)
        
        # cut to observation time
        Nt = len(hp)
        tm = np.arange(Nt) * self.dt
        T_sec = self.T * SIDEREALYEAR_J2000DAY*24*3600
        hp, hc, tm  = hp[tm<=T_sec], hc[tm<=T_sec], tm[tm<=T_sec]
        
        # interpolate if needed
        if approx_t:
            self.t, self.hp, self.hc = tm, hp, hc
        else:
            self._interp_hphc(tm, hp, hc, t, kind="spline")
        return (self.hp, self.hc)
        

if __name__ == "__main__":
    import doctest
    pEMRI_FSEF = dict({
        "MassMBHB": 1e6*un.Msun,
        "MassSOBHB": 10*un.Msun,
        "InitialSemiLatusRect": 12,
        "InitialEccentricity": 0.4,
        'InitialPolarAngleL':  3.9*un.rad,
        'InitialAzimuthalAngleL': 2.3535*un.rad,
        'InitialPhasePhiR': 0.123*un.rad,
        'InitialPhase': 0.456*un.rad,
        "Distance": 400*un.Mpc,
        "Cadence": 5.*un.s,
        "ObservationDuration": 3155814.97635456*un.s,
        "eps": 1e-5,
        "EclipticLatitude": 1.7*un.rad,
        "EclipticLongitude": 1.0471975511965976*un.rad,
    })
    doctest.testmod()
    
