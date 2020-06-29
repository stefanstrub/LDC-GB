import doctest
import numpy as np
import importlib
import tempfile

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

pMBHB = dict({'EclipticLatitude': 0.312414, #"radian"),
              'EclipticLongitude': -2.75291,# "radian"),
              'CoalescenceTime': 28086000.0,# 's'),
              'Distance':  9.14450149011798,# 'Gpc'),
              'InitialAzimuthalAngleL': 3.9,# 'radian'),
              'InitialPolarAngleL': 2.3535, #'radian'),
              'Mass1': 132628.202,# "SolarMass"),
              'Mass2': 30997.2481,# "SolarMass"),
              'PhaseAtCoalescence':  3.8, #'Radian'),
              'PolarAngleOfSpin1': 0.0,#'Radian'),
              'PolarAngleOfSpin2': 0.0,#'Radian'),
              'Redshift': 1.2687,# 'dimensionless'),
              'Spin1': 0.9481998052314212, #'MassSquared'),
              'Spin2': 0.9871324769575264, #'MassSquared'),
              'Cadence': 5.,
              'ObservationDuration':95.})#, 's') })

pEMRI = dict({"AzimuthalAngleOfSpin": 3.3, #'Radian'),
              "Distance": 4.066032714225506, #'Gpc'),
              "EclipticLatitude": -1.0865079083647613, #"Radian"),
              "EclipticLongitude": 5.45915569996,# 'Radian'),
              "InitialAlphaAngle": 5.144987211440423,# "Radian"),
              "InitialAzimuthalOrbitalFrequency": 0.0006865, #"Hz"),
              "InitialAzimuthalOrbitalPhase": 1.774141989722132,# "Radian"),
              "InitialEccentricity": 0.3365614512812651,# 'Unitless'),
              "InitialTildeGamma": 4.783451533651018,# 'Radian'),
              "LambdaAngle": 1.0976000000002033, #'Radian'),
              "MassOfCompactObject": 25.89699960277051,# 'SolarMass'),
              "MassOfSMBH": 1330605.9163026307,# 'SolarMass'),
              "PlungeTime": 41327816.5734, #'s'),
              "PolarAngleOfSpin": 2.3754411347141318, #'Radian'),
              "SMBHspin": 0.9671,# "MassSquared"),
              #'Cadence': (5., 's'),
})

tmp_filename = next(tempfile._get_candidate_names()) #pylint: disable=W0212


for mod in ["ldc.waveform.waveform",
            "ldc.waveform.source",
            "ldc.lisa.orbits",  "ldc.lisa.projection", "ldc.lisa.noise" , 
            "ldc.io.hdf5"]:
    module = importlib.import_module(mod)
    for sub in module.__all__:
        submodule = importlib.import_module(module.__name__+"."+sub)
        doctest.testmod(submodule, extraglobs=globals(), verbose=True)

