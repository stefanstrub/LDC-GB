import ldc.waveform.fastGB as FB
import numpy as np
#from memory_profiler import profile

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

def fourier2timedomain(XYZ, del_t, Tobs):
    Xf,Yf,Zf = XYZ
    Xt = np.fft.irfft(Xf)*(1.0/del_t)
    Yt = np.fft.irfft(Yf)*(1.0/del_t)
    Zt = np.fft.irfft(Zf)*(1.0/del_t)
    tim = np.arange(len(Xt))*del_t
    if (tim[-1] > Tobs): # cut the data at observation time
        i_end = np.argwhere(tim > Tobs)[0][0]
        Xt = Xt[:i_end]
        Yt = Yt[:i_end]
        Zt = Zt[:i_end]
        tim = tim[:i_end]
    return tim, Xt, Yt, Zt


del_t = 15
Tobs = 365*24*60*60*2
GB = FB.FastGB(delta_t=del_t, T=Tobs) # in seconds
freqT, X, Y, Z = GB.get_fd_tdixyz(template=pGB,
                                  oversample=4,
                                  simulator='synthlisa')

trange, X, Y, Z = fourier2timedomain((X,Y,Z), del_t, Tobs)

import FastGB as FB_MLDC
GB = FB_MLDC.FastGB("Test", dt=del_t, Tobs=Tobs, orbit="analytic")
bet = pGB["EclipticLatitude"]
lam = pGB["EclipticLongitude"]
Amp = pGB["Amplitude"]
f0 = pGB["Frequency"]
fdot = pGB["FrequencyDerivative"]
iota = pGB["Inclination"]
psi = pGB["Polarization"]
phi0 = pGB["InitialPhase"]
prm = np.array([f0, fdot, bet, lam, Amp, iota, psi, phi0])
X1, Y1, Z1 = GB.onefourier(simulator='synthlisa', params=prm, buffer=None,
                           T=Tobs, dt=del_t, algorithm='Michele', oversample=4)



trange1, X1, Y1, Z1 = fourier2timedomain((X1,Y1,Z1), del_t, Tobs)

import matplotlib.pyplot as plt
plt.figure()
plt.plot(trange1, X1)
plt.figure()
plt.plot(trange, X, alpha=0.5)
