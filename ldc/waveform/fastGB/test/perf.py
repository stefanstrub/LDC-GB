import numpy as np
import matplotlib.pyplot as plt

import ldc.waveform.fastGB as fastGB
from ldc.lisa.orbits import Orbits

import time

dt = 15
t_max = 60*60*24*365
orbits = Orbits.type(dict({'orbit_type':'analytic', 'nominal_arm_length':2.5e9,
                           "initial_position": 0, "initial_rotation": 0}))


pGB = dict({'Amplitude': 1.07345e-22,# "strain"
            'EclipticLatitude': 0.312414, # "radian"
            'EclipticLongitude': -2.75291,# "radian"
            'Frequency': 0.00135962, #"Hz"
            'FrequencyDerivative': 8.94581279e-19,# "Hz^2"
            'Inclination': 0.523599,# "radian"
            'InitialPhase': 3.0581565, #"radian"
            'Polarization': 3.5621656}) #"radian"

N = 100
t_c = np.zeros((N))
t_p = np.zeros((N))

GB = fastGB.FastGB(delta_t=dt, T=t_max, orbits=orbits) # in seconds
for i in range(N):
    t0 = time.time()
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator='synthlisa', radler=True)
    t1 = time.time()
    t_c[i] = t1-t0
    
for i in range(N):
    t0 = time.time()
    Xs2, Ys2, Zs2 = GB.get_fd_tdixyz(template=pGB, oversample=4, radler=False)
    t1 = time.time()
    t_p[i] = t1-t0

print(f"C version took {np.min(t_c)} s")
print(f"python version took {np.min(t_p)} s (x {np.round(np.min(t_p)/np.min(t_c),2)} slower)")


plt.figure()
i = 1
for tdi1, tdi2 in zip([Xs, Ys, Zs], [Xs2, Ys2, Zs2]):
    plt.subplot(3,1,i)
    plt.plot(tdi1.f, np.abs(tdi1), label="C")
    plt.plot(tdi2.f, np.abs(tdi2), label="python")
    i += 1
    plt.legend()
    plt.axis([pGB['Frequency']-1e-6, pGB['Frequency']+1e-6, None, None])
