import ldc.waveform.fastGB as FB
import numpy as np

pGB = dict({'Amplitude': 1.07345e-22,# "strain" 
            'EclipticLatitude': 0.312414, # "radian"
            'EclipticLongitude': -2.75291,# "radian" 
            'Frequency': 0.00135962, #"Hz"
            'FrequencyDerivative': 8.94581279e-19,# "Hz^2" 
            'Inclination': 0.523599,# "radian" 
            'InitialPhase': 3.0581565, #"radian"
            'Polarization': 3.5621656}) #"radian"
 
GB = FB.FastGB(delta_t=15, T=365*24*60*60) # in seconds
freqT, X, Y, Z = GB.get_fd_tdixyz(template=pGB,
                                  oversample=4,
                                  simulator='synthlisa')

