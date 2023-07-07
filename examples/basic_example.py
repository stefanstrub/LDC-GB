import numpy as np
import time

import unittest

try:
    import cupy as xp

    gpu_available = True

except (ImportError, ModuleNotFoundError) as e:
    import numpy as xp

    gpu_available = False

from gbgpu.gbgpu import GBGPU

from gbgpu.utils.constants import *

import subprocess as sp
import os


import matplotlib.pyplot as plt

from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries, window
import ldc.waveform.fastGB as fastGB

def get_gpu_memory():
    command = "nvidia-smi --query-gpu=memory.free --format=csv"
    memory_free_info = sp.check_output(command.split()).decode('ascii').split('\n')[:-1][1:]
    memory_free_values = [int(x.split()[0]) for i, x in enumerate(memory_free_info)]
    return memory_free_values

noise_model = "SciRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))

dt = 15.0
Tobs = 4.0 * YEAR
gb = GBGPU(use_gpu=gpu_available)

N = None
num_bin = 1*10**4
amp = 1e-22  # amplitude
f0 = 2e-3  # f0
fdot = 1e-16  # fdot
fddot = 0.0
phi0 = 0.1  # initial phase
iota = 0.2  # inclination
psi = 0.3  # polarization angle
lam = 0.4  # ecliptic longitude
beta_sky = 0.5  # ecliptic latitude
e1 = 0.2  # eccentricity of inner binary
beta1 = 0.5  # TODO: fill in
A2 = 19.5  # third body amplitude parameter
omegabar = 0.0  # omegabar parameter
e2 = 0.3  # eccentricity of third body
P2 = 0.6  # period of third body
T2 = 0.0  # time of periapsis passage of third body


# parameters
amp =  4.1204141965278374e-23  # amplitude
f0 = 0.004433557542242183  # f0
fdot = 7.39441302760258e-16 # fdot
fddot = 0.0 # fddot
phi0 = 2.295212550  # initial phase
iota = 0.96958009  # inclination
psi = 2.39042002 # polarization angle
lam = -1.5329 # ecliptic longitude
beta_sky =  0.1499696035197409  # ecliptic latitude

pGB = {'Amplitude': amp, 'Frequency': f0, 'FrequencyDerivative': fdot, 'InitialPhase': phi0, 'Inclination': iota, 'Polarization': psi, 'EclipticLongitude': lam, 'EclipticLatitude': beta_sky}

pGB = {'Amplitude': 1.1825788845002637e-21, 'EclipticLatitude': 1.1739903067531205, 'EclipticLongitude': -0.6355263029565035, 'Frequency': 0.0004948448022481871, 'FrequencyDerivative': -2.3640456997176845e-20, 'Inclination': 1.0413192157582787, 'InitialPhase': 3.8336721324711696, 'Polarization': 1.6675784744164501, 'IntrinsicSNR': 16.16742088451012}
# pGB = {'Amplitude': 3.563101227159202e-22, 'EclipticLatitude': 0.9967797425930816, 'EclipticLongitude': 0.6717384477032504, 'Frequency': 0.0018053741125663857, 'FrequencyDerivative': -6.465604466159056e-18, 'Inclination': 0.75113460407633, 'InitialPhase': 0.5908632463124097, 'Polarization': 1.1001953056658853, 'IntrinsicSNR': 103.32915328679616}
amp = pGB['Amplitude']
f0 = pGB['Frequency']
fdot = pGB['FrequencyDerivative']
phi0 = pGB['InitialPhase']
iota = pGB['Inclination']
psi = pGB['Polarization']
lam = pGB['EclipticLongitude']
beta_sky = pGB['EclipticLatitude']



amp_in = np.full(num_bin, amp)
f0_in = np.full(num_bin, f0)
fdot_in = np.full(num_bin, fdot)
fddot_in = np.full(num_bin, fddot)
phi0_in = np.full(num_bin, phi0)
iota_in = np.full(num_bin, iota)
psi_in = np.full(num_bin, psi)
lam_in = np.full(num_bin, lam)
beta_sky_in = np.full(num_bin, beta_sky)
e1_in = np.full(num_bin, e1)
beta1_in = np.full(num_bin, beta1)
A2_in = np.full(num_bin, A2)
P2_in = np.full(num_bin, P2)
omegabar_in = np.full(num_bin, omegabar)
e2_in = np.full(num_bin, e2)
T2_in = np.full(num_bin, T2)

modes = np.array([1, 2, 3])

length = int(Tobs / dt)

freqs = np.fft.rfftfreq(length, dt)
data_stream_length = len(freqs)

data = [
    1e-24 * xp.ones(data_stream_length, dtype=np.complex128),
    1e-24 * xp.ones(data_stream_length, dtype=np.complex128),
]

Nmodel = get_noise_model(noise_model, freqs)
Sn = Nmodel.psd(freq=freqs, option="X")
SA = Nmodel.psd(freq=freqs, option="A")
SE = Nmodel.psd(freq=freqs, option="E")
ST = Nmodel.psd(freq=freqs, option="T")

noise_factor = [
    xp.array(SA),
    xp.array(SA),
]

# noise_factor = [
#     1e-28 * xp.ones(data_stream_length, dtype=np.float64),
#     1e-28 * xp.ones(data_stream_length, dtype=np.float64),
# ]

params_circ = np.array(
    [amp_in, f0_in, fdot_in, fddot_in, phi0_in, iota_in, psi_in, lam_in, beta_sky_in,]
)
params_circ[1,0] = 2.0001e-3
params_circ[1,1] = 2.001e-3
params_circ[1,2] = 2.01e-3
params_circ[1,3] = 2.1e-3
params_circ[1,3] = 3e-3

num = 10

#######################
####  CIRCULAR ########
#######################
A_inj, E_inj = gb.inject_signal(
    amp,
    f0,
    fdot,
    fddot,
    phi0,
    iota,
    psi,
    lam,
    beta_sky,
    N=N,
    dt=dt,
    T=Tobs,
)

data = [xp.array(A_inj), xp.array(E_inj)]

gb.d_d = 0
st = time.perf_counter()
for _ in range(num):
    print(get_gpu_memory())
    like = gb.get_ll(
        params_circ, data, noise_factor, N=N, dt=dt, T=Tobs,
    )
et = time.perf_counter()
print("circ:", (et - st) / num, "per binary:", (et - st) / (num * num_bin))


df = 7e-6
GB = fastGB.FastGB(delta_t=dt, T=Tobs)
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator="synthlisa")
Af = (Zs - Xs)/np.sqrt(2.0)
Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
indexes = np.logical_and(freqs >= Af.f.values[0], freqs <= Af.f.values[-1]) 



SNR2 = np.sum( np.real(data[0][indexes].get() * np.conjugate(Af.data) + data[1][indexes].get() * np.conjugate(Ef.data))/noise_factor[0][indexes].get() )
hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /noise_factor[0][indexes].get())

logliks = 4.0*Xs.df*( SNR2 - 0.5 * hh )

plt.figure()
plt.semilogy(freqs,np.abs(data[0].get()))
plt.xlim(f0 - df, f0 + df)
plt.show()

plt.figure()
plt.semilogy(freqs,np.abs(data[0].get()))
plt.semilogy(freqs,np.abs(A_inj),'.')
plt.semilogy(Af.f,np.abs(Af))
plt.semilogy(Af.f,np.abs(data[0][indexes].get()))
plt.xlim(f0 - df, f0 + df)
plt.show()

et = time.perf_counter()
print("circ:", (et - st) / num, "per binary:", (et - st) / (num * num_bin))


#######################
####  ECCENTRIC #######
#######################


params_ecc = np.array(
    [
        amp_in,
        f0_in,
        fdot_in,
        fddot_in,
        phi0_in,
        iota_in,
        psi_in,
        lam_in,
        beta_sky_in,
        e1_in,
        beta1_in,
    ]
)

modes = np.array([1, 2, 3, 4])
A_inj, E_inj = gb.inject_signal(
    amp,
    f0,
    fdot,
    fddot,
    phi0,
    iota,
    psi,
    lam,
    beta_sky,
    e1,
    beta1,
    modes=modes,
    N=N,
    dt=dt,
    T=Tobs,
)

st = time.perf_counter()
for _ in range(num):
    like = gb.get_ll(params_ecc, data, noise_factor, N=N, dt=dt, modes=modes, T=Tobs,)
et = time.perf_counter()
print(
    "ecc ({} modes):".format(len(modes)),
    (et - st) / num,
    "per binary:",
    (et - st) / (num * num_bin),
)


##############################
## CIRCULAR / THIRD BODY #####
##############################

params_circ_third = np.array(
    [
        amp_in,
        f0_in,
        fdot_in,
        fddot_in,
        phi0_in,
        iota_in,
        psi_in,
        lam_in,
        beta_sky_in,
        A2_in,
        omegabar_in,
        e2_in,
        P2_in,
        T2_in,
    ]
)

A_inj, E_inj = gb.inject_signal(
    amp,
    f0,
    fdot,
    fddot,
    phi0,
    iota,
    psi,
    lam,
    beta_sky,
    A2,
    omegabar,
    e2,
    P2,
    T2,
    modes=np.array([2]),
    N=N,
    dt=dt,
    T=Tobs,
)

st = time.perf_counter()
for _ in range(num):
    like = gb.get_ll(
        params_circ_third, data, noise_factor, N=N, dt=dt, modes=np.array([2]), T=Tobs,
    )
et = time.perf_counter()
print("circ / third:", (et - st) / num, "per binary:", (et - st) / (num * num_bin))

##############################
## ECCENTRIC / THIRD BODY #####
##############################

params_full = np.array(
    [
        amp_in,
        f0_in,
        fdot_in,
        fddot_in,
        phi0_in,
        iota_in,
        psi_in,
        lam_in,
        beta_sky_in,
        e1_in,
        beta1_in,
        A2_in,
        omegabar_in,
        e2_in,
        P2_in,
        T2_in,
    ]
)

modes = np.array([1, 2, 3, 4])
A_inj, E_inj = gb.inject_signal(
    amp,
    f0,
    fdot,
    fddot,
    phi0,
    iota,
    psi,
    lam,
    beta_sky,
    e1,
    beta1,
    A2,
    omegabar,
    e2,
    P2,
    T2,
    modes=modes,
    N=N,
    dt=dt,
    T=Tobs,
)

st = time.perf_counter()
for _ in range(num):
    like = gb.get_ll(params_full, data, noise_factor, N=N, dt=dt, modes=modes, T=Tobs,)
et = time.perf_counter()
print(
    "ecc/third ({} modes):".format(len(modes)),
    (et - st) / num,
    "per binary:",
    (et - st) / (num * num_bin),
)

breakpoint()
