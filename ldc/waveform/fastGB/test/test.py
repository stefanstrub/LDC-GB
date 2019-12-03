import ldc.waveform.fastGB as FB
import numpy as np

pGB = dict({'Amplitude': 1.07345e-22, 
            'EclipticLatitude': 0.312414, 
            'EclipticLongitude': -2.75291,
            'Frequency': 0.00135962, 
            'FrequencyDerivative': 8.94581279e-19, 
            'Inclination': 0.523599 ,
            'InitialPhase': 3.0581565,
            'Polarization': 3.5621656})
 
GB = FB.FastGB(template=pGB, T=24*60*60)
x, y, z = GB.get_fd_tdixyz()
print(x,y,z)
