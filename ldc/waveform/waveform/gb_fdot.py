""" Compute waveforms h+ and hx for galactic binaries in time domain. """


import numpy as np
from astropy import units as un
from ldc.waveform.waveform.hphc import HpHc


#pylint:disable=E1101

class GB_fdot(HpHc):
    """Compute waveforms h+ and hx of a galactic binary

    Vectorized sources are supported in this case, by setting vectors
    for each parameter.
    """

    parameter_map = {'amplitude': 'Amplitude',
                     'phi0': 'InitialPhase',
                     'f': lambda p: np.array([p['Frequency']]),
                     'dfdt': lambda p: np.array([p['FrequencyDerivative']])}

    def precomputation(self):
        """ Additional precomputation. """
        super().precomputation()
        self.cos_inc = np.cos(self.source_parameters['Inclination'])

    def display(self):
        """ Display the source and precomputed parameters. """
        super().display()
        print("Internal parameters:")
        print('- cos(inc)  =', self.cos_inc, 'rad')

    def info(self):
        """ Return GB parameter names and units
        """
        units = {'EclipticLatitude':         'rad',
                 'EclipticLongitude':        'rad',
                 'Amplitude':                '1',
                 'Frequency':                'Hz',
                 'FrequencyDerivative':      'Hz2',
                 'Inclination':              'rad',
                 'Polarization':             'rad',
                 'InitialPhase':             'rad'}
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

        >>> GB = HpHc.type("my-galactic-binary", "GB", "TD_fdot")
        >>> hp,hc = GB.compute_hphc_td(np.arange(0,100,10), pGB)
        >>> print(hp[0:3], hc[0:3] )
        [1.13239259e-22 1.00151885e-22 8.63340609e-23] [1.49869582e-22 1.58865610e-22 1.66702965e-22]
        """
        if source_parameters is not None:
            self.set_param(source_parameters)

        # Check the approximant and call appropriate function
        if self.approximant == 'TD_fdot':
            phase = -self.phi0 + np.pi*(2.*t[:, None]*self.f[None, :] \
                                        + (t*t)[:, None]*self.dfdt[None, :])
            phase = phase.squeeze()
            hSp = -np.cos(phase)*self.amplitude * (1 + self.cos_inc * self.cos_inc)
            hSc = -np.sin(phase)*2.*self.amplitude * self.cos_inc
        else:
            raise NotImplementedError

        hp, hc = self.source2SSB(hSp, hSc) # Convert to SSB
        if not set_attr:
            return hp, hc
        else:
            self.hp, self.hc = hp, hc
            self.t = t
            return self.hp, self.hc


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

    doctest.testmod()
