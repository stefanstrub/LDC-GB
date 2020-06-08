import doctest
import numpy as np
import importlib

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

for mod in ["ldc.waveform.waveform"]:
    module = importlib.import_module(mod)
    for sub in module.__all__:
        submodule = importlib.import_module(module.__name__+"."+sub)
        doctest.testmod(submodule, extraglobs=globals(), verbose=True)

