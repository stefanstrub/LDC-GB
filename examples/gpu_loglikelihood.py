# import torch
# from torch.utils.data import Dataset, DataLoader, random_split
# from torch import nn
# from torch import distributions

import matplotlib.pyplot as plt

import numpy as np
import time

# Making code agnostic to CPU/GPU
def std_get_wrapper(arg):
    return arg

def cuda_get_wrapper(arg):
    return arg.get()

try:
   import cupy as cp
   gpu = True
   get_wrapper = cuda_get_wrapper
except:
   import numpy as cp
   gpu = False
   get_wrapper = std_get_wrapper

print('GPU',gpu)

from gbgpu.gbgpu import GBGPU
from gbgpu.noisemodel import AnalyticNoise 

# Size of the training set.
num_samples = 10**4

# We choose a very narrow frequency range.
f0_lower = 0.010062
f0_upper = 0.010084

# Amplitude range.
amp_lower = 1e-23
amp_upper = 1e-21

# Sample f0 and amp from a uniform prior.
f0 = cp.random.uniform(f0_lower, f0_upper, num_samples)
amp = cp.random.uniform(amp_lower, amp_upper, num_samples)

# Fixed parameters
ones = cp.ones(num_samples)
fdot = 1.79e-15 * ones
lam  = 4.36 * ones
beta = 2.18 * ones
iota = 0.67 * ones
phi0 = 5.48 * ones
psi  = 0.43 * ones

# Package parameters into arrays.

sampling_parameters = cp.vstack((f0, amp)).T
all_parameters = cp.vstack((amp, f0, fdot, cp.zeros(num_samples), -phi0, iota, psi, lam, beta)).T

# Initialise waveform generator.
gb = GBGPU(use_gpu=gpu)

# Waveform settings
Tobs = 31536000.0*4  # One-year observation
dt = 15.0  # Sample rate (Nyquist is safely larger than the maximum frequency we will encounter)
df = 1./Tobs
N_points = 128

# Generate the waveforms.
start = time.time()
gb.run_wave(*all_parameters.T, N = N_points, dt = dt, T = Tobs, oversample = 1)
print('wave generation time', time.time()-start)
start = time.time()
gb.run_wave(*all_parameters.T, N = N_points, dt = dt, T = Tobs, oversample = 1)
print('wave generation time', time.time()-start)

f_min = 0.010059
f_max = 0.0100861

# Define the frequency grid.
num_bins = int((f_max - f_min) / df) + 1
sample_frequencies = cp.linspace(f_min, f_max, num=num_bins)

noise = AnalyticNoise(sample_frequencies)
psd_A, psd_E = noise.psd(option="A"), noise.psd(option="E")


asd_A = cp.sqrt(psd_A)
asd_E = cp.sqrt(psd_E)
psd = cp.zeros((2, len(psd_A)))
psd[0] = psd_A
psd[1] = psd_E

k_min = round(f_min/df)
k_max = round(f_max/df)
num = len(sample_frequencies)


# These indices describe how to stitch the waveform into the larger frequency grid.
i_start = (get_wrapper(gb.start_inds) - k_min).astype(cp.int32)
i_end = (get_wrapper(gb.start_inds) - k_min + gb.N).astype(cp.int32)

# PyTorch by default uses float32, and that should be sufficient for our purposes.
# Here we use complex64 since the frequency-domain strain is complex.

# start = time.time()
# A_whitened = cp.zeros((num_samples, num-len(gb.A[0])), dtype=cp.complex128)
# # E_whitened = cp.zeros((num_samples, num-len(gb.A[0])), dtype=cp.complex128)
# A_whitened_ext = cp.concatenate((gb.A,A_whitened), axis=1)

# start = time.time()
# for i in range(num_samples):
#     A_whitened_ext[i] = cp.roll(A_whitened_ext[i],i_start[i])
#     # E_whitened_ext[i] = cp.roll(E_whitened_ext[i],i_start[i], axis=1)
# print('time',time.time()-start)

# start = time.time()
# A_whitened_ext_2 = [cp.roll(A_whitened_ext[i],i) for i in i_start]
# print('time',time.time()-start)

# start = time.time()
# A_whitened = cp.zeros((num_samples, num), dtype=np.complex128)
# E_whitened = cp.zeros((num_samples, num), dtype=np.complex128)
# for i in range(num_samples):
#     A_whitened[i,i_start[i]:i_end[i]] = gb.A[i]
#     E_whitened[i,i_start[i]:i_end[i]] = gb.E[i]
# print('arangement time cuda',time.time()-start)

# start = time.time()
# A_whitened *= cp.sqrt(4 * df) / asd_A
# E_whitened *= cp.sqrt(4 * df) / asd_E
# print('noise time',time.time()-start)

start = time.time()
numpy_A = gb.A.get()
numpy_E = gb.E.get()
A_whitened = np.zeros((num_samples, num), dtype=np.complex128)
E_whitened = np.zeros((num_samples, num), dtype=np.complex128)
for i in range(num_samples):
    A_whitened[i,i_start[i]:i_end[i]] = numpy_A[i]
    E_whitened[i,i_start[i]:i_end[i]] = numpy_E[i]
A_whitened = cp.asarray(A_whitened)
E_whitened = cp.asarray(E_whitened)
print('arangement time numpy',time.time()-start)

noise = cp.random.normal(size = A_whitened[0].real.shape).astype(cp.float32)
dataA = cp.zeros(num, dtype=cp.complex128)
dataA.real = A_whitened[0].real + noise[0]
noise = cp.random.normal(size = A_whitened[0].real.shape).astype(cp.float32)
dataA.imag = A_whitened[0].imag + noise[0]

noise = cp.random.normal(size = E_whitened[0].real.shape).astype(cp.float32)
dataE = cp.zeros(num, dtype=cp.complex128)
dataE.real = E_whitened[0].real + noise[0]
noise = cp.random.normal(size = E_whitened[0].real.shape).astype(cp.float32)
dataE.imag = E_whitened[0].imag + noise[0]


start = time.time()
SNR2 = np.sum( np.real(dataA * np.conjugate(A_whitened) + dataE * np.conjugate(E_whitened))/asd_A, axis=1 )
hh = np.sum((np.absolute(A_whitened)**2 + np.absolute(E_whitened)**2) /asd_A, axis=1 )
loglikelihood = 4.0*df*( SNR2 - 0.5 * hh )
print('loglikelihood time',time.time()-start)

start = time.time()
gb.get_ll(*all_parameters.T, data=)
print('noise time',time.time()-start)


start = time.time()
A_whitened *= cp.sqrt(4 * df) / asd_A
E_whitened *= cp.sqrt(4 * df) / asd_E
print('noise time',time.time()-start)


# Plot a sample waveform
plt.plot(get_wrapper(sample_frequencies), get_wrapper(A_whitened[0].real))
# plt.plot(get_wrapper(sample_frequencies), get_wrapper(dataA.real))
plt.xscale('log')
plt.xlabel('$f$')
plt.ylabel('Re $A$')
plt.show()

# For best training, parameters should be standardized (zero mean, unit variance across the training set).

parameters_mean = np.mean(sampling_parameters, axis=0)
parameters_std = np.std(sampling_parameters, axis=0)

parameters_standardized = (sampling_parameters - parameters_mean) / parameters_std
parameters_standardized = parameters_standardized.astype(np.float32)

# Prepare the (complex) frequency-domain data for the (real) neural network. To do so we simply concatenate the real and imaginary parts into an array of doubled length.

waveforms = cp.hstack((A_whitened.real, A_whitened.imag))

class WaveformDataset(Dataset):
    
    def __init__(self, parameters, waveforms):
        self.parameters = parameters
        self.waveforms = waveforms

    def __len__(self):
        return len(self.parameters)

    def __getitem__(self, idx):
        params = self.parameters[idx]
        signal = self.waveforms[idx]
        
        # Add unit normal noise to the signal
        noise = cp.random.normal(size = signal.shape).astype(cp.float32)
        data = signal + noise
        
        return torch.tensor(data), torch.tensor(params)

waveform_dataset = WaveformDataset(parameters_standardized, waveforms)

# We can sample from the WaveformDataset. This gives us pairs of data and parameters, different noise realizations each time.

x, y = waveform_dataset[2]
plt.plot(x)
plt.show()


# start = time.time()
# zeros_left = []
# for i in range(num_samples):
#     zeros_left.append(cp.zeros( i_start[i], dtype=cp.complex128))
# zeros_left = np.asarray(zeros_left)
# A_whitened_ext = cp.concatenate((zeros_left,gb.A), axis=1)


# A_whitened_ext = cp.concatenate([cp.roll(A_whitened_ext[i],i) for i in i_start])
#     E_whitened[i,i_start[i]:i_end[i]] = gb.E[i]

data = cp.zeros((2, len(A_whitened[0])))
data[0] = A_whitened[0]
data[1] = E_whitened[0]


# gb.d_d = 
# start = time.time()
# gb.get_ll(all_parameters.T, data=data, psd = psd, N = N_points, dt = dt, T = Tobs, oversample = 1)
# print('time',time.time()-start)

start = time.time()
A_whitened_slow = cp.empty((num_samples, num), dtype=cp.complex64)
E_whitened_slow = cp.empty((num_samples, num), dtype=cp.complex64)

for i in range(num_samples):
    x = cp.zeros(num, dtype=cp.complex128)
    x[i_start[i]:i_end[i]] = gb.A[i]
    x *= cp.sqrt(4 * df) / asd_A
    A_whitened_slow[i] = x
    
    x = cp.zeros(num, dtype=cp.complex128)
    x[i_start[i]:i_end[i]] = gb.E[i]
    x *= cp.sqrt(4 * df) / asd_E
    E_whitened_slow[i] = x
print('time',time.time()-start)
print('end')