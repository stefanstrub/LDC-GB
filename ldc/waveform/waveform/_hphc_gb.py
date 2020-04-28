import numpy as np

def info(self):
    """ Return GB parameter names and units
    """
    units = {'EclipticLatitude':         'Radian',
             'EclipticLongitude':        'Radian',
             'Amplitude':                'strain',
             'Frequency':                'Hz',
             'FrequencyDerivative':      'Hz^2',
             'Inclination':              'Radian',
             'Polarization':             'Radian',
             'InitialPhase':             'Radian'}
    return units


def check_param(self):
    """ Check parameters and their units
    """
    for k in list(self.info().keys()):
        assert k in self.pnames
    assert self.units["Polarization"].lower() in ["radian", "rad", "r"]
    assert self.units["EclipticLatitude"].lower() in ["radian", "rad", "r"]
    assert self.units["EclipticLongitude"].lower() in ["radian", "rad", "r"]
    assert self.units["InitialPhase"].lower() in ["radian", "rad", "r"]
    assert self.units["Inclination"].lower() in ["radian", "rad", "r"]
    assert self.units["Frequency"].lower() in ["hertz", "hz"]

def compute_hphc_td(self, t, source_parameters=None, approx_t=False, set_attr=False):
    """ Return hp,hx for a time samples in t.

    Source parameters can be updated at the same time. 
    Keyword approx_t is ignore for here. 

    >>> GB = HpHcGB("my-galactic-binary", "GB", "TD_fdot")
    >>> hp,hc = GB.compute_hphc_td(np.arange(0,100,10), pGB)
    >>> print(hp[0:3], hc[0:3] )
    [1.13239259e-22 1.00151885e-22 8.63340609e-23] [1.49869582e-22 1.58865610e-22 1.66702965e-22]
    """
    if source_parameters != None:
        self.set_param(source_parameters)

    # Check the approximant and call appropriate function
    if self.approximant == 'TD_fdot':
        phase = -self.phi0 + np.pi*(2.*t[:, None]*self.f[None, :] \
                                    + (t*t)[:, None]*self.dfdt[None, :])
        phase = phase.squeeze()
        hSp = -np.cos(phase)*self.amplitude * (1 + self.cosInc * self.cosInc)
        hSc = -np.sin(phase)*2.*self.amplitude * self.cosInc
    else:
        raise NotImplementedError

    hp, hc = self.source2SSB(hSp, hSc) # Convert to SSB
    if not set_attr:
        return hp,hc
    else:
        self.hp, self.hc = hp, hc
        self.t = t
        return self.hp, self.hc


if __name__ == "__main__":
    import doctest
    import numpy as np
    from hphc import HpHcGB

    pGB = dict({'Amplitude': 1.07345e-22,#, "strain"), 
                'EclipticLatitude': 0.312414,#, "radian"),
                'EclipticLongitude': -2.75291,# "radian"), 
                'Frequency': 0.00135962,# "Hz"),
                'FrequencyDerivative': 8.94581279e-19,# "Hz^2"), 
                'Inclination': 0.523599 ,# "radian"), 
                'InitialPhase': 3.0581565,# "radian"), 
                'Polarization': 3.5621656})#,# "radian")})

    
    doctest.testmod()
