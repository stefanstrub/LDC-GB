import doctest
import numpy as np
import importlib
import tempfile
from astropy import units as un
from ldc.waveform.waveform import HpHc
from ldc.lisa.orbits import Orbits
import lisaorbits

my_orbits = lisaorbits.EqualArmlengthOrbits(dt=100, size=100)
my_orbits.write(f'orbits.h5')

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

pMBHB = dict({'Mass1':                  132628.202,
                'Mass2':                  30997.2481,
                'Spin1':                  0.9481998052314212, 
                'Spin2':                  0.9871324769575264, 
                'EclipticLatitude':       0.312414,
                'EclipticLongitude':      -2.75291,
                'Redshift':               2.0,
                'Distance':               15974.456786495544,
                'Cadence':                5.,
                'ObservationDuration':    95.,
                'CoalescenceTime':        28086000.0,
                'InitialAzimuthalAngleL': 3.9,
                'InitialPolarAngleL':     2.3535,
                'PhaseAtCoalescence':     3.8,
                'PolarAngleOfSpin1':      0.0,
                'PolarAngleOfSpin2':      0.0,
                })
pSOBHB = dict({"Mass1":                 50.,
                "Spin1":                 0.0,
                "Mass2":                 40.,
                "Spin2":                 0.0,
                "EclipticLatitude":      1.7,
                "EclipticLongitude":     1.0471975511965976,
                "Inclination":           1.0471975511965976,
                "InitialFrequency":      1.2e-2,
                "InitialPhase":          0.7,
                "Polarization":          1.2,
                "Redshift":              2.0,
                "Distance":              15974.456786495544,
                'Cadence':               5,                   
                'ObservationDuration':   189348898.58127362/6})

pEMRI = dict({"AzimuthalAngleOfSpin": 3.3*un.rad,
              "Distance": 4.066032714225506*un.Gpc,
              "EclipticLatitude": -1.0865079083647613*un.rad,
              "EclipticLongitude": 5.45915569996*un.rad,
              "InitialAlphaAngle": 5.144987211440423*un.rad,
              "InitialAzimuthalOrbitalFrequency": 0.0006865*un.Hz,
              "InitialAzimuthalOrbitalPhase": 1.774141989722132*un.rad,
              "InitialEccentricity": 0.3365614512812651,
              "InitialTildeGamma": 4.783451533651018*un.rad,
              "LambdaAngle": 1.0976000000002033*un.rad,
              "MassOfCompactObject": 25.89699960277051*un.Msun,
              "MassOfSMBH": 1330605.9163026307*un.Msun,
              "PlungeTime": 41327816.5734*un.s,
              "PolarAngleOfSpin": 2.3754411347141318*un.rad,
              "SMBHspin": 0.9671,# "MassSquared"),
})

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

config = {'nominal_arm_length':2.5e9, #"m", 
          'initial_rotation':0, #'rad', 
          'initial_position':0,#'rad')]
          'orbit_type':'analytic'}
GB = HpHc.type("my-galactic-binary", "GB", "TD_fdot")
GB.set_param(pGB)
orbits = Orbits.type(config)

tmp_filename = next(tempfile._get_candidate_names()) #pylint: disable=W0212

if __name__ == "__main__":
    for mod in ["ldc.waveform.waveform",
                "ldc.waveform.source",
                "ldc.lisa.orbits",  "ldc.lisa.projection", "ldc.lisa.noise" , 
                "ldc.io.hdf5", "ldc.common.series"]:
        module = importlib.import_module(mod)
        for sub in module.__all__:
            submodule = importlib.import_module(module.__name__+"."+sub)
            doctest.testmod(submodule, extraglobs=globals(), verbose=True,
                            raise_on_error=True)

