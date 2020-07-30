import doctest
import numpy as np
import importlib
import tempfile
from astropy import units as un
from ldc.waveform.waveform import HpHc
from ldc.lisa.orbits import Orbits

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

pSOBBH = dict({"AzimuthalAngleOfSpin1": 0.0*un.rad,
               "AzimuthalAngleOfSpin2": 0.0*un.rad,
               "Distance": 0.8217407069275701*un.Gpc,
               "EclipticLatitude": 0.23339632679489664*un.rad,
               "EclipticLongitude": 1.1798*un.rad,
               "Inclination": 1.1508*un.rad,
               "InitialFrequency": 0.0074076*un.Hz,
               "InitialPhase": 1.2622*un.rad,
               "Mass1": 31.033*un.Msun,
               "Mass2": 19.918*un.Msun,
               "PolarAngleOfSpin1": 2.7329*un.rad,
               "PolarAngleOfSpin2": 2.2947*un.rad,
               "Polarization": 3.7217*un.rad,
               "Redshift": 0.16454,
               "Spin1": 0.4684,# "MassSquared"),
               "Spin2": 0.979,# "MassSquared"),
               'Cadence': 5.*un.s,
               'ObservationDuration': 95.*un.s})

pMBHB = dict({'EclipticLatitude': 0.312414*un.rad,
              'EclipticLongitude': -2.75291*un.rad,
              'CoalescenceTime': 28086000.0*un.s,
              'Distance':  9.14450149011798*un.Gpc,
              'InitialAzimuthalAngleL': 3.9*un.rad,
              'InitialPolarAngleL': 2.3535*un.rad,
              'Mass1': 132628.202*un.Msun,
              'Mass2': 30997.2481*un.Msun,
              'PhaseAtCoalescence':  3.8*un.rad,
              'PolarAngleOfSpin1': 0.0*un.rad,
              'PolarAngleOfSpin2': 0.0*un.rad,
              'Redshift': 1.2687,
              'Spin1': 0.9481998052314212, #'MassSquared'),
              'Spin2': 0.9871324769575264, #'MassSquared'),
              'Cadence': 5.*un.s,
              'ObservationDuration':95.*un.s})

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

config = {'nominal_arm_length':2.5e9, #"m", 
          'initial_rotation':0, #'rad', 
          'initial_position':0,#'rad')]
          'orbit_type':'analytic'}

GB = HpHc.type("my-galactic-binary", "GB", "TD_fdot")
GB.set_param(pGB)
orbits = Orbits.type(config)

tmp_filename = next(tempfile._get_candidate_names()) #pylint: disable=W0212


for mod in ["ldc.waveform.waveform",
            "ldc.waveform.source",
            "ldc.lisa.orbits",  "ldc.lisa.projection", "ldc.lisa.noise" , 
            "ldc.io.hdf5"]:
    module = importlib.import_module(mod)
    for sub in module.__all__:
        submodule = importlib.import_module(module.__name__+"."+sub)
        doctest.testmod(submodule, extraglobs=globals(), verbose=True)

