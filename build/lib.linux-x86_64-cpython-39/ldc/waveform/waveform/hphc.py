""" Compute waveforms h+ and hx for all types of sources at a given time. """

from abc import ABC, abstractmethod
import copy
import importlib

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from astropy import units as un

try:
    import cupy as cp
except ImportError:
    cp = None

#pylint:disable=C0103
#pylint:disable=W0201

def interp_cuda_dev(x_vals, x, y):
    """
    Interpolate [x, y] and evaluate in x_vals
    Input and output are cupy arrays  
    """
    return cp.interp(x_vals, x, y)
    
def _interp(t, x, tnew, kind='spline', der=False, integr=False):
    """ Perform the interpolation of hx, hp at time samples in t. """
    if cp is not None and cp.get_array_module(tnew) is cp:
        # if tnew is in device, use cuda interpolation
        # t, x are in host, tnew is already in device (see method interp_hphc)
        return interp_cuda_dev(tnew, cp.array(t), cp.array(x))
    if kind == 'spline':
        if der:
            splh = spline(t, x).derivative()
        elif integr:
            splh = spline(t, x).antiderivative()
        else:
            splh = spline(t, x)
        return (splh(tnew), splh)
    #hp = np.interp(tnew, t, x)
    itp = interp1d(t, x, kind=kind)
    return itp(tnew), itp


class HpHc(ABC):
    """Compute waveforms h+ and hx for all types of sources.

    An HpHc is identified by a:

    - source type in [MBHB, EMRI, GB, SOBBH, SGWB]
    - approximant in [IMRPhenomD, AK, PhenomD, None]
    - set of parameters.

    The set of parameters depends on the approximant and the source type.
    """

    parameter_map = {}

    def __getattr__(self, parname):
        if parname in self.parameter_map:
            parmapping = self.parameter_map[parname]
            return (self.source_parameters[parmapping] if isinstance(parmapping, str)
                    else parmapping(self.source_parameters))
        elif parname in self.__dict__.keys():
            # this will throw the right exception if we don't have this attribute
            return self.__dict__[parname] #if isinstance(parname)
        raise AttributeError

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
            module = importlib.import_module("ldc.waveform.waveform.bhb_imrphenomD")
            return getattr(module, "BHB_IMRPhenomD")(source_name, source_type, approximant)
        if source_type == "SBBH" and approximant == "IMRPhenomD":
            module = importlib.import_module("ldc.waveform.waveform.bhb_imrphenomD")
            return getattr(module, "BHB_IMRPhenomD")(source_name, source_type, approximant)
        if source_type == "GB":
            module = importlib.import_module("ldc.waveform.waveform.gb_fdot")
            return getattr(module, "GB_fdot")(source_name, source_type, approximant)

        # SOBHB models
        if source_type == "SOBHB" and approximant == "IMRPhenomD":
            module = importlib.import_module("ldc.waveform.waveform.bhb_imrphenomD")
            # return getattr(module, "SOBHB_IMRPhenomD")(source_name, source_type, approximant)
            return getattr(module, "BHB_IMRPhenomD")(source_name, source_type, approximant)

        # EMRI models
        if source_type == "EMRI" and approximant == "AK":
            module = importlib.import_module("ldc.waveform.waveform.emri_ak")
            return getattr(module, "EMRI_AK")(source_name, source_type, approximant)
        if source_type == "EMRI" and approximant == "FSEF":
            module = importlib.import_module("ldc.waveform.waveform.emri_fsef")
            return getattr(module, "EMRI_FastSchwarschildEccentricFlux")(source_name, source_type, approximant)
        raise ValueError("Invalid source_type %s (approximant=%s)"%(source_type, approximant))

    @property
    def pnames(self):
        """ Shortcut to parameter name """
        if isinstance(self.source_parameters, dict):
            return list(self.source_parameters.keys())
        return list(self.source_parameters.dtype.names)

    def split(self):
        """Return a list of HpHc object, one for each set of parameter in
        current object.

        >>> hphc = HpHc.type("demo", "GB", "None")
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
        for i, GW in enumerate(GWS):
            pars = [(k, p[k][i]) if isinstance(p[k], np.ndarray) else (k, p[k])
                    for k in self.pnames]
            GW.source_parameters = copy.copy(self.source_parameters)
            GW.set_param(dict(pars))
        return GWS

    def set_param(self, param, units=None):
        """Set or update all source parameters

        Args:
            param (dict, array, record array): if simple array, alphabetical order is assumed
            units (dict): units can also be given through param using astropy quantities.

        >>> HpHc = HpHc.type("demo", "GB", "None")
        >>> HpHc.set_param(pGB)
        >>> HpHc.display()
        Source parameters:
        Amplitude : 1.07345e-22  [ 1 ]
        EclipticLatitude : 0.312414  [ rad ]
        EclipticLongitude : -2.75291  [ rad ]
        Frequency : 0.00135962  [ Hz ]
        FrequencyDerivative : 8.94581279e-19  [ Hz2 ]
        Inclination : 0.523599  [ rad ]
        InitialPhase : 3.0581565  [ rad ]
        Polarization : 3.5621656  [ rad ]
        Internal parameters:
        - cos(inc)  = 0.8660252915835662 rad

        """
        if isinstance(param, dict):
            self.source_parameters = param
        elif (isinstance(param, np.recarray) or (isinstance(param, np.record)) or
              (isinstance(param, np.ndarray) and
               hasattr(param.dtype, 'names'))):
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

    def set_units(self, units=None):
        """ Set and convert to default units.
        """
        default_units = self.info()
        for k in self.pnames:
            if isinstance(self.source_parameters[k], un.Quantity):
                conv = self.source_parameters[k].to(un.Unit(default_units[k]))
                self.source_parameters[k] = conv.value
            elif units is not None and k in units and un.Unit(units[k])!=un.dimensionless_unscaled:
                qty = self.source_parameters[k]*un.Unit(units[k])
                qty = qty.to(un.Unit(default_units[k]))
                self.source_parameters[k] = qty.value
        self.units = default_units


    @abstractmethod
    def check_param(self):
        """ Check parameters and their units
        """
        pass

    @abstractmethod
    def compute_hphc_td(self, t, source_parameters=None, approx_t=False, set_attr=False):
        """Return hp,hx for a time samples in t and store t,hp,hx as
        attributes.

        Source parameters can be updated at the same time.  if
        approx_t is True, no interpolation is done, and actual time
        range used might differ from the one given.

        """
        pass

    def interp_hphc(self, t, cuda=False):
        """ Interpolate hp,hx on a given time range.

        A first call to compute_hphc_td is assumed here.
        """
        if cuda is True:
            if cp.get_array_module(t) is np:
                # Move array to device if it is not
                cu_t = cp.array(t)
            else:
                cu_t = t
            hp = _interp(self.t, self.hp, cu_t)
            hc = _interp(self.t, self.hc, cu_t) 
            return hp, hc
        if "i_hp" in self.__dict__:
            return self.i_hp(t), self.i_hc(t)
        hp, self.i_hp = _interp(self.t, self.hp, t)
        hc, self.i_hc = _interp(self.t, self.hc, t)
        return hp, hc

    def _interp_hphc(self, tm, hp, hc, t, **kwargs):
        """ Interpolate hp, hc and set them as attr.

        Extrapolation is a 0 padding.
        """
        t_start, t_end = max(tm[0], t[0]), min(tm[-1], t[-1])
        i_st = np.argwhere(t >= t_start)[0][0]
        i_en = np.argwhere(t >= t_end)[0][0]
        self.hp, self.hc = np.zeros(len(t)), np.zeros(len(t))
        self.hp[i_st:i_en+1], self.i_hp = _interp(tm, hp, t[i_st:i_en+1], **kwargs)
        self.hc[i_st:i_en+1], self.i_hc = _interp(tm, hc, t[i_st:i_en+1], **kwargs)
        self.t = t

    def precomputation(self):
        """ Precompute projectors
        """
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
        if isinstance(sin_d, np.ndarray) and sin_d.ndim > 0:
            u = np.array([sin_l, -cos_l, np.zeros((len(sin_l)))])
        else:
            u = np.array([sin_l, -cos_l, 0])
        self.basis = k, v, u

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

class NumericHpHc(HpHc):
    """ h+ and hx from arrays.

    """

    def __init__(self, t, hp, hc, lon, lat, name="num-hphc"):
        """ Set original values and interpolant
        """
        self.t = t
        self.hp = hp
        self.hc = hc
        self.i_hp = spline(t, hp)
        self.i_hc = spline(t, hc)
        self.source_type = "Numeric"
        self.source_name = name
        self.set_param(dict({"EclipticLatitude":lat,
                             "EclipticLongitude":lon}))

    def compute_hphc_td(self, t):
        return interp_hphc(t)

    def check_param(self):
        """ Check parameters and their units
        """
        for k in list(self.info().keys()):
            assert k in self.pnames
        assert self.units["EclipticLatitude"].lower() in ["radian", "rad", "r"]
        assert self.units["EclipticLongitude"].lower() in ["radian", "rad", "r"]

    def info(self):
        """ Return GB parameter names and units
        """
        units = {'EclipticLatitude':         'rad',
                 'EclipticLongitude':        'rad'}
        return units

if __name__ == "__main__":
    import doctest

    pGB = dict({'Amplitude': 1.07345e-22,
                'EclipticLatitude': 0.312414*un.rad,
                'EclipticLongitude': -2.75291*un.rad,
                'Frequency': 0.00135962*un.Hz,
                'FrequencyDerivative': 8.94581279e-19*un.Unit('Hz2'),
                'Inclination': 0.523599*un.rad,
                'InitialPhase': 3.0581565*un.rad,
                'Polarization': 3.5621656*un.rad})

    pGBl = dict({'Amplitude': np.array([1.07345e-22, 1.0125e-22]),
                 'EclipticLatitude': np.array([0.312414, 1.015463])*un.rad,
                 'EclipticLongitude': np.array([-2.75291, 0.512364])*un.rad,
                 'Frequency': np.array([0.00135962, 0.0005478])*un.Hz,
                 'FrequencyDerivative': np.array([8.94581279e-19,
                                                  8.45279e-19]*un.Unit('Hz2')),
                 'Inclination': np.array([0.523599, 0.15548])*un.rad,
                 'InitialPhase': np.array([3.0581565, 3.546841])*un.rad,
                 'Polarization': np.array([3.5621656, 3.124485])*un.rad})


    doctest.testmod()
