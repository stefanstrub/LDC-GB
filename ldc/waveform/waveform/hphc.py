""" Compute waveforms h+ and hx for all types of sources at a given time. """

from abc import ABC, abstractmethod
import copy
import importlib

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as spline

#pylint:disable=C0103
#pylint:disable=W0201

def _interp(t, x, tnew, kind='spline', der=False, integr=False):
    """ Perform the interpolation of hx, hp at time samples in t. """
    if kind == 'spline':
        if der:
            splh = spline(t, x).derivative()
        elif integr:
            splh = spline(t, x).antiderivative()
        else:
            splh = spline(t, x)
        return (splh(tnew), splh)
    else:
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
        else:
            # this will throw the right exception if we don't have this attribute
            return self.__dict__[parname] #if isinstance(parname)


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
            module = importlib.import_module("ldc.waveform.waveform.bbh_imrphenomD")
            return getattr(module, "BBH_IMRPhenomD")(source_name, source_type, approximant)
        elif source_type == "GB":
            module = importlib.import_module("ldc.waveform.waveform.gb_fdot")
            return getattr(module, "GB_fdot")(source_name, source_type, approximant)
        elif source_type == "EMRI":
            module = importlib.import_module("ldc.waveform.waveform.emri_ak")
            return getattr(module, "EMRI_AK")(source_name, source_type, approximant)
        elif source_type == "SOBBH":
            module = importlib.import_module("ldc.waveform.waveform.sobbh_phenomD")
            return getattr(module, "SOBBH_phenomD")(source_name, source_type, approximant)
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
        if isinstance(self.source_parameters, dict):
            return list(self.source_parameters.keys())
        else:
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

    # def add_param(self, param, value, units='dimensionless'):
    #     """ add or update parameter
    #     """
    #     self.source_parameters[param] = value
    #     self.units[param] = units


    def set_param(self, param, units='default'):
        """Set or update all source parameters

        param can be
        - a dict
        - a vector of values (alphabetical order is assumed)
        - a record array

        >>> HpHc = HpHc.type("demo", "GB", "None")
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
            self.source_parameters = param
        elif (isinstance(param, np.recarray) or
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

    def set_units(self, units='default'):
        """ Set units to default value.
        """
        if units != "default":
            raise ValueError("Units can't be changed for now")
        self.units = self.info()

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

    def interp_hphc(self, t):
        """ Interpolate hp,hx on a given time range.

        A first call to compute_hphc_td is assumed here.
        """
        if "i_hp" in self.__dict__:
            return self.i_hp(t), self.i_hc(t)
        else:
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
        if isinstance(sin_d, np.ndarray):
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


if __name__ == "__main__":
    import doctest

    pGB = dict({'Amplitude': 1.07345e-22,#, "strain"),
                'EclipticLatitude': 0.312414,#, "radian"),
                'EclipticLongitude': -2.75291,# "radian"),
                'Frequency': 0.00135962,# "Hz"),
                'FrequencyDerivative': 8.94581279e-19,# "Hz^2"),
                'Inclination': 0.523599,# "radian"),
                'InitialPhase': 3.0581565,# "radian"),
                'Polarization': 3.5621656})#,# "radian")})

    pGBl = dict({'Amplitude': np.array([1.07345e-22, 1.0125e-22]), #"strain"),
                 'EclipticLatitude': np.array([0.312414, 1.015463]),# "radian"),
                 'EclipticLongitude': np.array([-2.75291, 0.512364]),# "radian"),
                 'Frequency': np.array([0.00135962, 0.0005478]), #"Hz"),
                 'FrequencyDerivative': np.array([8.94581279e-19, 8.45279e-19]), #"Hz^2"),
                 'Inclination': np.array([0.523599, 0.15548]), #"radian"),
                 'InitialPhase': np.array([3.0581565, 3.546841]),# "radian"),
                 'Polarization': np.array([3.5621656, 3.124485])})#, "Radian"),

    doctest.testmod()
