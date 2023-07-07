import numpy as np
import time

import matplotlib.pyplot as plt

from gbgpu.gbgpu import GBGPU
from gbgpu.thirdbody import GBGPUThirdBody

from gbgpu.utils.constants import *
from gbgpu.utils.utility import *

gb = GBGPU(use_gpu=True)

dt = 10.0
Tobs = 4.0 * YEAR


# number of points in waveform
# if None, will determine inside the code based on amp, f0 (and P2 if running third-body waveform)
N = None

# number of binaries to batch
num_bin = 10

# parameters
amp = 2e-23  # amplitude
f0 = 2e-3  # f0
fdot = 7.538331e-18  # fdot
fddot = 0.0 # fddot
phi0 = 0.1  # initial phase
iota = 0.2  # inclination
psi = 0.3  # polarization angle
lam = 0.4  # ecliptic longitude
beta_sky = 0.5  # ecliptic latitude


# parameters
amp =  4.1204141965278374e-23  # amplitude
# f0 = 0.004433557542242183  # f0
fdot = 7.39441302760258e-16 # fdot
fddot = 0.0 # fddot
phi0 = 2.295212550  # initial phase
iota = 0.96958009  # inclination
psi = 2.39042002 # polarization angle
lam = -1.5329 # ecliptic longitude
beta_sky =  0.1499696035197409  # ecliptic latitude

# for batching
amp_in = np.full(num_bin, amp)
f0_in = np.full(num_bin, f0)
fdot_in = np.full(num_bin, fdot)
fddot_in = np.full(num_bin, fddot)
phi0_in = np.full(num_bin, phi0)
iota_in = np.full(num_bin, iota)
psi_in = np.full(num_bin, psi)
lam_in = np.full(num_bin, lam)
beta_sky_in = np.full(num_bin, beta_sky)

pGB = {'Amplitude': amp, 'Frequency': f0, 'FrequencyDerivative': fdot, 'InitialPhase': phi0, 'Inclination': iota, 'Polarization': psi, 'EclipticLongitude': lam, 'EclipticLatitude': beta_sky}

params = np.array(
    [amp_in, f0_in, fdot_in, fddot_in, phi0_in, iota_in, psi_in, lam_in, beta_sky_in,]
)

gb.run_wave(*params, N=N, dt=dt, T=Tobs, oversample=2)

from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries, window
import ldc.waveform.fastGB as fastGB
GB = fastGB.FastGB(delta_t=dt, T=Tobs)
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator="synthlisa")
Af = (Zs - Xs)/np.sqrt(2.0)

# signal from first binary
A = gb.A[0]
X = gb.X[0]
freqs = gb.freqs[0]
print("signal length:", A.shape)
plt.plot(freqs.get(), A.get())
plt.plot(Af.f, Af)
plt.ylabel("TDI-A Channel", fontsize=16)
plt.xlabel("freqs (Hz)", fontsize=16)
dx = 7e-7
plt.xlim(f0 - dx, f0 + dx)
plt.show()

X = gb.X[0]
freqs = gb.freqs[0]
print("signal length:", A.shape)
plt.plot(freqs.get(), X.get())
plt.plot(Xs.f, Xs)
plt.ylabel("TDI-X Channel", fontsize=16)
plt.xlabel("freqs (Hz)", fontsize=16)
dx = 7e-7
plt.xlim(f0 - dx, f0 + dx)
plt.show()



gb.get_ll(*params, )

breakpoint()