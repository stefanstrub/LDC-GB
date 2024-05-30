import ldc.waveform.fastGB as fastGB

from ldc.lisa import orbits
from astropy import units as un
import numpy as np
import matplotlib.pyplot as plt

pGB = dict({'Amplitude': 1.07345e-22,# "strain"
            'EclipticLatitude': 0.312414, # "radian"
            'EclipticLongitude': -2.75291,# "radian"
            'Frequency': 0.00135962, #"Hz"
            'FrequencyDerivative': 8.94581279e-19,# "Hz^2"
            'Inclination': 0.523599,# "radian"
            'InitialPhase': 3.0581565, #"radian"
            'Polarization': 3.5621656}) #"radian"

del_t = 15
Tobs = 365*24*60*60*0.5
lisa_orbits = orbits.Orbits.type(dict({"nominal_arm_length":2.5e6*un.km,
                                       "initial_rotation":0*un.rad,
                                       "initial_position":0*un.rad,
                                       "orbit_type":"analytic"}))
trange = np.arange(0, Tobs, del_t)
N = fastGB.get_buffersize(Tobs, pGB["Frequency"], pGB["Amplitude"], oversample=4)

GB_C = fastGB.FastGB(delta_t=del_t, T=Tobs, orbits=lisa_orbits) #historical / radler version
GB = fastGB.VFastGB(delta_t=del_t, T=Tobs, orbits=lisa_orbits, N=N) 

# in chronological order
X1t, Y1t, Z1t = GB_C.get_td_tdixyz(template=pGB, oversample=4, radler=True)
X2t, Y2t, Z2t = GB_C.get_td_tdixyz(template=pGB, oversample=4, new=False)
X3t, Y3t, Z3t = GB_C.get_td_tdixyz(template=pGB, oversample=4)
X4t, Y4t, Z4t = GB.get_td_tdixyz(template=pGB)

plt.figure()
for j,(T1,T2,T3,T4) in enumerate([(X1t, X2t, X3t, X4t),
                                  (Y1t, Y2t, Y3t, Y4t),
                                  (Z1t, Z2t, Z3t, Z4t)]):
    plt.subplot(3, 2, 2*j+1)
    T1.plot(label="LDC radler")
    plt.plot(trange, T2, label="LDC sangria", alpha=0.5)
    plt.legend(loc="lower right")
    plt.subplot(3, 2, 2*j+2)
    plt.plot(trange, np.array(T1)-T2, label="difference")
    plt.plot(trange, T2-T3, label="cross check 1")
    plt.plot(trange, T3-T4, label="cross check 2")
    plt.legend()

# in chronological order
X1f, Y1f, Z1f = GB_C.get_fd_tdixyz(template=pGB, oversample=4, radler=True)
X2f, Y2f, Z2f = GB_C.get_fd_tdixyz(template=pGB, oversample=4, new=False)
X3f, Y3f, Z3f = GB_C.get_fd_tdixyz(template=pGB, oversample=4)
X4f, Y4f, Z4f = GB.get_fd_tdixyz(template=pGB)

plt.figure()
for j,(T1,T2,T3,T4) in enumerate([(X1t, X2t, X3t, X4t),
                                  (Y1t, Y2t, Y3t, Y4t),
                                  (Z1t, Z2t, Z3t, Z4t)]):
    plt.subplot(3, 2, 2*j+1)
    plt.plot(np.abs(T1), label="LDC radler")
    plt.plot(np.abs(T2), label="LDC sangria", alpha=0.5)
    plt.legend(loc="lower right")
    plt.subplot(3, 2, 2*j+2)
    plt.plot(np.abs(T1-T2), label="difference")
    plt.plot(np.abs(T2-T3), label="cross check 1")
    plt.plot(np.abs(T3-T4), label="cross check 2")
    plt.legend()
    
number = 1000
import timeit
print(timeit.timeit('GB_C.get_fd_tdixyz(template=pGB, oversample=4, radler=True, xarray=False)',
                    globals=globals(), number=number)/number)
print(timeit.timeit('GB_C.get_fd_tdixyz(template=pGB, oversample=4, new=False, xarray=False)',
                    globals=globals(), number=number)/number)
print(timeit.timeit('GB_C.get_fd_tdixyz(template=pGB, oversample=4, xarray=False)',
                    globals=globals(), number=number)/number)
print(timeit.timeit('GB.get_fd_tdixyz(template=pGB, xarray=False)',
                    globals=globals(), number=number)/number)



[f0, fdot, ampl, theta, phi, psi, incl, phi0] = GB_C._parse_template(pGB)

X3f, Y3f, Z3f = GB_C.get_fd_tdixyz(f0=f0, fdot=fdot, ampl=ampl, theta=theta, phi=phi,
                                   psi=psi, incl=incl, phi0=phi0, oversample=4)
X4f, Y4f, Z4f = GB.get_fd_tdixyz(f0=f0, fdot=fdot, ampl=ampl, theta=theta, phi=phi,
                                 psi=psi, incl=incl, phi0=phi0)

plt.figure()
for j,(T3,T4) in enumerate([(X3f, X4f),
                            (Y3f, Y4f),
                            (Z3f, Z4f)]):
    plt.subplot(3, 2, 2*j+1)
    plt.plot(T3.f, np.abs(T3), label="LDC")
    plt.plot(T4.f, np.abs(T4), label="wrapper", alpha=0.5)
    plt.legend(loc="lower right")
    plt.subplot(3, 2, 2*j+2)
    plt.plot(np.abs(T3-T4), label="cross check")
