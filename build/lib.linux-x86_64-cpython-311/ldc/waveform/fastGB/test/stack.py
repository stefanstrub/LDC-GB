import numpy as np
import matplotlib.pyplot as plt
import ldc.waveform.fastGB as FB
from ldc.common.series import XYZ2AET, TDI, FrequencySeries

#default template
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
GB = FB.FastGB(delta_t=del_t, T=Tobs)

N = 10 # number of sources
Nf = int(((Tobs/del_t)-1)/2)
Xs = FrequencySeries(np.zeros((Nf), dtype=np.complex128), df=1/Tobs, kmin=0)
Ys = FrequencySeries(np.zeros((Nf), dtype=np.complex128), df=1/Tobs, kmin=0)
Zs = FrequencySeries(np.zeros((Nf), dtype=np.complex128), df=1/Tobs, kmin=0)
for i in range(N):

    pGBi = pGB.copy()
    pGBi["Frequency"] += np.random.randn()*5e-6

    X,Y,Z = GB.get_fd_tdixyz(template=pGBi, oversample=4)
    kmin = X.attrs['kmin']
    Xs[kmin:kmin+len(X)] += X 
    Ys[kmin:kmin+len(X)] += Y 
    Zs[kmin:kmin+len(X)] += Z


# Look in freq domain around 13mHz
plt.figure()
plt.subplot(211)
selec = Xs.sel(f=slice(0.0012, 0.0014))
plt.plot(selec.f, np.abs(selec))

# Go to time domain
plt.subplot(212)
Xt = Xs.ts.ifft(dt=del_t)
plt.plot(Xt.t, Xt)

    
