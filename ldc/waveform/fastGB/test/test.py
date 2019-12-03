import ldc.waveform.fastGB as FB
from LISAhdf5 import LISAhdf5, ParsUnits
import numpy as np

pGB = dict({'Amplitude': (1.07345e-22, "strain"), 
            'EclipticLatitude': (0.312414, "radian"),
            'EclipticLongitude': (-2.75291, "radian"), 
            'Frequency': (0.00135962, "Hz"),
            'FrequencyDerivative': (8.94581279e-19, "Hz^2"), 
            'Inclination': (0.523599 , "radian"), 
            'InitialPhase': (3.0581565, "radian"), 
            'Polarization': (3.5621656, "radian"),
            'Cadence': (15., 's')})
 
pGB = ParsUnits(dict([(k,v[0]) for k,v in pGB.items()]),
                dict([(k,v[1]) for k,v in pGB.items()]))
pars = pGB.pars

GB = FB.FastGB(template=pars, T=24*60*60)
x, y, z = GB.get_fd_tdixyz()
print(x,y,z)
