""" Compute waveforms h+ and hx for all types of sources at a given time. """

from abc import ABC, abstractmethod
import copy

import numpy as np
import pyfftw
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as spline

from ldc.common import constants
from ldc.common import tools
#pylint:disable=C0103

class HpHc(ABC):
    """Compute waveforms h+ and hx for all types of sources.

    An HpHc is identified by a:
    - source type in [MBHB, EMRI, GB, SOBBH, SGWB]
    - approximant in [IMRPhenomD, AK, PhenomD, None]
    - set of parameters.

    The set of parameters depends on the approximant and the source type.
    """
    def __init__(self, source_name, source_type, approximant):
        """ Initialization common to all sources. """
        self.source_name = source_name
        self.source_type = source_type
        self.approximant = approximant
        self.reference_frame = "SSB" # default
        
    @classmethod
    def type(cls, source_name, source_type, approximant):
        """ Return instance corresponding to source type. """
        if source_type == "MBHB" or (source_type is None and approximant == "IMRPhenomD"):
            return HpHcMBHB(source_name, source_type, approximant)
        elif source_type == "GB":
            return HpHcGB(source_name, source_type, approximant)
        elif source_type == "EMRI":
            return HpHcEMRI(source_name, source_type, approximant)
        elif source_type == "SOBBH":
            return HpHcSOBBH(source_name, source_type, approximant)
        raise ValueError("Invalid source_type %s (approximant=%s)"%(source_type, approximant))

        
    @staticmethod
    def from_file(filename, source_name=None, index=None):
        """Return the hxhp instance corresponding to source type read from
        file, or passed as argument.

        If source_name is None, take the first source.
        hx, hx and time samples are read from file.

        TODO: remove dependancy to MLDC
        """
        from LISAhdf5 import LISAhdf5
        h5 = LISAhdf5(filename, mode='r')
        if source_name is None and index is not None:
            source_name = h5.getSourcesName()[index]
        param = h5.getSourceParameters(source_name)
        source_type = param.get("SourceType")
        approximant = param.get("Approximant")

        hphc = HpHc.type(source_name, source_type, approximant)
        hphc.set_param(param.pars, units=param.units)
        try:
            hphc.hp, hphc.hx, hphc.t = h5.getSourceHpHc(source_name)
        except:
            pass
        return hphc

    @property
    def pnames(self):
        """ Shortcut to parameter name """
        return list(self.source_parameters.keys()) if isinstance(self.source_parameters, dict) else list(self.source_parameters.dtype.names)
    
    def split(self):
        """Return a list of HpHc object, one for each set of parameter in
        current object.
        
        >>> hphc = HpHcGB("demo", "GB", "None")
        >>> hphc.set_param(pGBl)
        >>> print(hphc.eclLat)
        [0.312414 1.015463]
        >>> hphc_list = hphc.split()
        >>> print(hphc_list[0].eclLat)
        0.312414
        """
        p = self.source_parameters
        ns = [len(p[k]) if isinstance(p[k], np.ndarray) else 1
              for k in self.pnames]
        ns = np.array(ns).max()
        GWS = [HpHc.type(self.source_name, self.source_type, self.approximant) for s in range(ns)]
        for i,GW in enumerate(GWS):
            pars = [(k,p[k][i]) if isinstance(p[k], np.ndarray) else (k,p[k])
                    for k in self.pnames]
            GW.source_parameters = copy.copy(self.source_parameters)
            GW.set_param(dict(pars))
        return GWS

    def add_param(self, param, value, units='dimensionless'):
        """ add or update parameter 
        """
        self.source_parameters[param] = value
        self.units[param] = units
        
    
    def set_param(self, param, units='default'):
        """Set or update all source parameters

        param can be 
        - a dict 
        - a vector of values (alphabetical order is assumed)
        - a record array

        >>> HpHc = HpHcGB("demo", "GB", "None")
        >>> HpHc.set_param(pGB)
        >>> HpHc.display()
        Source parameters:
        Amplitude : 1.07345e-22  [ strain ]
        EclipticLatitude : 0.312414  [ Radian ]
        EclipticLongitude : -2.75291  [ Radian ]
        Frequency : 0.00135962  [ Hz ]
        FrequencyDerivative : 8.94581279e-19  [ Hz^2 ]
        Inclination : 0.523599  [ Radian ]
        InitialPhase : 3.0581565  [ Radian ]
        Polarization : 3.5621656  [ Radian ]
        Internal parameters:
        - phi0 =  3.0581565 rad
        - f    =  [0.00135962] Hz
        - dfdt =  [8.94581279e-19] Hz/s
        - amplitude =  3.0581565
        - cos(inc)  = 0.8660252915835662
        """
        if isinstance(param, dict):
            self.source_parameters = param #np.rec.fromarrays(param.values(),names=list(param.keys()))
        elif isinstance(param, np.recarray) or (isinstance(param, np.ndarray) and hasattr(param.dtype, 'names')):
            self.source_parameters = param
        elif hasattr(param, "__len__"):
            names = sorted(self.source_parameters.keys())
            self.source_parameters = dict(zip(names, param))
        else:
            raise ValueError("Invalid param")

        ## Pre-computation specific to the source
        self.set_units(units)
        self.check_param()
        self.precomputation()

    def set_units(self, units='default'):
        self.units = self.info() if units=='default' else units

    @abstractmethod
    def check_param(self):
        """ Check parameters and their units
        """ 
        pass
    
    @abstractmethod
    def compute_hphc_td(self, t, source_parameters=None, approx_t=False):
           """Return hp,hx for a time samples in t and store t,hp,hx as
           attributes.
           
           Source parameters can be updated at the same time.  if
           approx_t is True, no interpolation is done, and actual time
           range used might differ from the one given.

           """
           pass

    def interp_hphc(self, t, kind="spline"):
        """ Interpolate hp,hx on a given time range. 

        A first call to compute_hphc_td is assumed here. 
        """
        if hasattr(self, "i_hp"):
            return self.i_hp(t), self.i_hc(t)
        else:
            hp, self.i_hp = self._interp(self.t, self.hp, t)
            hc, self.i_hc = self._interp(self.t, self.hc, t)
            return hp, hc

    def _interp_hphc(self, tm, hp, hc, t, **kwargs):
        """ Interpolate hp, hc and set them as attr. 

        Extrapolation is a 0 padding. 
        """
        t_start, t_end = max(tm[0], t[0]), min(tm[-1], t[-1])
        i_st = np.argwhere(t>=t_start)[0][0]
        i_en = np.argwhere(t>= t_end)[0][0]
        self.hp, self.hc = np.zeros(len(t)), np.zeros(len(t))
        self.hp[i_st:i_en+1],self.i_hp = self._interp(tm, hp, t[i_st:i_en+1], **kwargs)
        self.hc[i_st:i_en+1],self.i_hc = self._interp(tm, hc, t[i_st:i_en+1], **kwargs)
        self.t = t
    
    def _interp(self, t, x, tnew, kind='spline', der=False, integr=False):
        """ Perform the interpolation of hx, hp at time samples in t. """
        if kind=='spline':
            if (der):
                splh =  spline(t, x).derivative()
            elif (integr):
                splh = splins(t, x).antiderivative()
            else:
                splh = spline(t, x)
            return (splh(tnew), splh)
        else:
            hp = np.interp(tnew, t, x)
            itp = interp1d(t, x, kind=kind)
            return itp(tnew), itp


    def precomputation(self):
        """ """
        
        p = self.source_parameters
        if 'Polarization' in self.pnames:
            self.pol = p['Polarization']
            self.cos2Pol = np.cos(2.*self.pol)
            self.sin2Pol = np.sin(2.*self.pol)
        self.eclLat = p['EclipticLatitude']
        self.eclLon = p['EclipticLongitude']

        ## Projection basis
        sin_d, cos_d = np.sin(self.eclLat), np.cos(self.eclLat)
        sin_l, cos_l = np.sin(self.eclLon), np.cos(self.eclLon)
        k = np.array([-cos_d * cos_l, -cos_d * sin_l, -sin_d])
        v = np.array([-sin_d * cos_l, -sin_d * sin_l, cos_d])
        u = np.array([sin_l, -cos_l, np.zeros((len(sin_l)))]) if isinstance(sin_d, np.ndarray) else np.array([sin_l, -cos_l, 0])
        self.basis = k,v,u

    def to_file(self, filename):
        """ Save hp, hx, t and source type and parameters to file. 

        TODO: remove dependancy to MLDC
        """
        from LISAhdf5 import LISAhdf5
        h5 = LISAhdf5(filename)
        from LISAhdf5 import ParsUnits
        units = self.source_parameters.copy()
        for k,v in self.units.items():
            units[k] = v
        pu = ParsUnits(pars_i=self.source_parameters, units_i=units)
        h5.addSource(self.source_name, pu,
                     overwrite=True, hphcData=np.vstack([self.t, self.hp, self.hc]).T)
            

    def source2SSB(self, hSp, hSc):
        """ Convert h+, hx from source frame to Solar System Barycenter.
            Convention defined in the LDC documentation. 
        """
        hp = hSp * self.cos2Pol - hSc * self.sin2Pol
        hc = hSp * self.sin2Pol + hSc * self.cos2Pol
        return hp, hc

    def display(self):
        """ Display the source parameters.
        """
        print("Source parameters:")
        for k in self.pnames:
            if k in self.units.keys():
                print(k, ":", self.source_parameters[k], " [", self.units[k], "]")
            else:
                print(k, ":", self.source_parameters[k])


class HpHcGB(HpHc):
    """Compute waveforms h+ and hx of a galactic binary
    
    Vectorized sources are supported in this case, by setting vectors
    for each parameter.
    """
    
    def precomputation(self):
        """ Additional precomputation. """
        super().precomputation()
        self.cosInc = np.cos(self.source_parameters['Inclination'])

    def display(self):
        """ Display the source and precomputed parameters. """
        super().display()
        print("Internal parameters:")
        print('- phi0 = ', self.phi0, 'rad')
        print('- f    = ', self.f, 'Hz')
        print('- dfdt = ', self.dfdt, 'Hz/s')
        print('- amplitude = ', self.phi0)
        print('- cos(inc)  =', self.cosInc)
        
    @property
    def amplitude(self):
        return self.source_parameters['Amplitude']
    @property
    def phi0(self):
        return self.source_parameters['InitialPhase']
    @property
    def f(self):
        return np.array([self.source_parameters['Frequency']])
    @property
    def dfdt(self):
        return np.array([self.source_parameters['FrequencyDerivative']])

    from ._hphc_gb import info, check_param, compute_hphc_td

class HpHcMBHB(HpHc):
    """ Compute waveforms h+ and hx of a massive black hole binary. 

    Vectorized sources are not supported in this case. 
    """

    @property
    def redshift(self):
        return self.source_parameters['Redshift']
    @property
    def phi0(self):
        return self.source_parameters['PhaseAtCoalescence']
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
    def m1s(self):
        return self.source_parameters['Mass1']
    @property
    def m2s(self):
        return self.source_parameters['Mass2']
    @property
    def tc(self):
        return self.source_parameters['CoalescenceTime']
    @property
    def dt(self):
        return self.source_parameters['Cadence']
    @property
    def Tobs(self):
        return self.source_parameters['ObservationDuration']
    
    def precomputation(self):
        """ Load required parameters and convert them in expected units. """
        super().precomputation()
        p = self.source_parameters
        theL = p['InitialPolarAngleL']
        phiL = p['InitialAzimuthalAngleL']
        self.pol, self.incl = tools.aziPolAngleL2PsiIncl(self.eclLat, self.eclLon, theL, phiL)
        self.cos2Pol = np.cos(2.*self.pol)
        self.sin2Pol = np.sin(2.*self.pol)
        self.a1 = np.cos(self.Stheta1s)*self.chi1s # For PhenomD we will use projections
        self.a2 = np.cos(self.Stheta2s)*self.chi2s
        self.dist = self.DL*1e6*constants.Nature.PARSEC_METER
        self.set_FD()
        
    from ._hphc_mbhb import info, check_param, set_FD, compute_hphc_td, IMRPhenomD_MBHB
    from ._hphc_mbhb import _IMRPhenomD_waveform
    
class HpHcSOBBH(HpHc):
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
        
    from ._hphc_sobbh import info, check_param, phenomD_SOBBH, compute_hphc_td, set_FD

class HpHcEMRI(HpHc):
    """ Compute waveforms h+ and hx of an EMRI.

    Vectorized sources are not supported in this case. 
    """

    def precomputation(self):
        """ Load required parameters and convert them in expected units. """
        import EMRI_AK
        super().precomputation()
        p = self.source_parameters
        #self.dt = p.getConvert('Cadence',LC.convT,'sec')
        self.wv = EMRI_AK.AK("zoom")
        self.wv.EclipticLatitude =  p['EclipticLatitude']
        self.wv.EclipticLongitude = p['EclipticLongitude']
        self.wv.PolarAngleOfSpin = p['PolarAngleOfSpin']
        self.wv.AzimuthalAngleOfSpin = p['AzimuthalAngleOfSpin']
        self.wv.Spin = p['SMBHspin']
        self.wv.MassOfCompactObject = p['MassOfCompactObject']
        self.wv.MassOfSMBH = p['MassOfSMBH']
        self.wv.InitialAzimuthalOrbitalFrequency = p['InitialAzimuthalOrbitalFrequency']
        self.wv.InitialAzimuthalOrbitalPhase = p['InitialAzimuthalOrbitalPhase']
        self.wv.InitialEccentricity = p['InitialEccentricity']
        self.wv.InitialTildeGamma = p['InitialTildeGamma']
        self.wv.InitialAlphaAngle = p['InitialAlphaAngle']
        self.wv.LambdaAngle = p['LambdaAngle']
        self.wv.Distance = p['Distance']*1.e9 # Gpc -> pc 

        
    from ._hphc_emri import info, check_param, compute_hphc_td

if __name__ == "__main__":
    import numpy as np
    import doctest

    pGB = dict({'Amplitude': 1.07345e-22,#, "strain"), 
                'EclipticLatitude': 0.312414,#, "radian"),
                'EclipticLongitude': -2.75291,# "radian"), 
                'Frequency': 0.00135962,# "Hz"),
                'FrequencyDerivative': 8.94581279e-19,# "Hz^2"), 
                'Inclination': 0.523599 ,# "radian"), 
                'InitialPhase': 3.0581565,# "radian"), 
                'Polarization': 3.5621656})#,# "radian")})
    
    pGBl = dict({'Amplitude': np.array([1.07345e-22, 1.0125e-22]), #"strain"), 
                 'EclipticLatitude': np.array([0.312414, 1.015463]),# "radian"),
                 'EclipticLongitude': np.array([-2.75291, 0.512364]),# "radian"), 
                 'Frequency': np.array([0.00135962, 0.0005478]), #"Hz"),
                 'FrequencyDerivative': np.array([8.94581279e-19, 8.45279e-19]), #"Hz^2"), 
                 'Inclination': np.array([0.523599, 0.15548 ]), #"radian"), 
                 'InitialPhase': np.array([3.0581565, 3.546841]),# "radian"), 
                 'Polarization': np.array([3.5621656, 3.124485])})#, "Radian"),

    doctest.testmod()
