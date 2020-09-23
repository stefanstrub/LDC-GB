import ldc.waveform.fastGB as FB
import numpy as np
#from memory_profiler import profile
import FastGB as FB_MLDC
import matplotlib.pyplot as plt


pGB = dict({'Amplitude': 1.07345e-22,# "strain" 
            'EclipticLatitude': 0.312414, # "radian"
            'EclipticLongitude': -2.75291,# "radian" 
            'Frequency': 0.00135962, #"Hz"
            'FrequencyDerivative': 8.94581279e-19,# "Hz^2" 
            'Inclination': 0.523599,# "radian" 
            'InitialPhase': 3.0581565, #"radian"
            'Polarization': 3.5621656}) #"radian"

#@profile
def loop(prm, Tobs, del_t):
    GB = FB.FastGB(delta_t=15, T=365*24*60*60) # in seconds
    
    for i in range(10000):
        freqT, X, Y, Z = GB.get_fd_tdixyz(template=pGB,
                                          oversample=4,
                                          simulator='synthlisa')
    return 0

#loop(prm, Tobs, del_t)

del_t = 15
Tobs = 365*24*60*60*0.5
GB = FB.FastGB(delta_t=del_t, T=Tobs) # in seconds

X, Y, Z = GB.get_td_tdixyz(template=pGB,
                           oversample=4,
                           simulator='synthlisa')
trange = np.arange(0, Tobs, del_t)

GB2 = FB.FastGB(delta_t=del_t, T=Tobs, ldc_orbits=1) # in seconds

X2, Y2, Z2 = GB2.get_td_tdixyz(template=pGB,
                               oversample=4,
                               simulator='synthlisa')

plt.figure()
for j,(T1,T2) in enumerate([(X, X2), (Y,Y2), (Z,Z2)]):
    plt.subplot(3, 2, 2*j+1)
    T1.plot(label="default orbits")
    T2.plot(label="LDC orbits", alpha=0.5)
    plt.legend(loc="lower right")
    plt.subplot(3, 2, 2*j+2)
    plt.plot(trange, np.array(T1)-np.array(T2), label="difference")


# GB = FB_MLDC.FastGB("Test", dt=del_t, Tobs=Tobs, orbit="analytic")
# bet = pGB["EclipticLatitude"]
# lam = pGB["EclipticLongitude"]
# Amp = pGB["Amplitude"]
# f0 = pGB["Frequency"]
# fdot = pGB["FrequencyDerivative"]
# iota = pGB["Inclination"]
# psi = pGB["Polarization"]
# phi0 = pGB["InitialPhase"]
# prm = np.array([f0, fdot, bet, lam, Amp, iota, psi, -phi0])
# X1, Y1, Z1 = GB.onefourier(simulator='synthlisa', params=prm, buffer=None,
#                            T=Tobs, dt=del_t, algorithm='Michele', oversample=4)

# X1t = X1.ifft(del_t) 
# Y1t = Y1.ifft(del_t)
# Z1t = Z1.ifft(del_t)

# plt.figure()
# for j,(T1,T2) in enumerate([(X, X1t), (Y,Y1t), (Z,Z1t)]):
#     plt.subplot(3, 2, 2*j+1)
#     T1.plot(label="newLDC")
#     plt.plot(trange, T2, label="MLDC", alpha=0.5)
#     plt.legend(loc="lower right")
#     plt.subplot(3, 2, 2*j+2)
#     plt.plot(trange, np.array(T1)-T2, label="difference")


