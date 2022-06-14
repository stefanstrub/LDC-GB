# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import scipy as sp

plt.rcParams["font.family"] = "serif"
plt.rcParams.update({'font.size': 15})

# Make two noise-free time series ---------------------------------
f_target = 0.0001
b_target = 0.0000001
a_target = 30
delta_target = 400000
period = 1 / f_target
sample_step_size = 15
sample_rate = 1 / sample_step_size
measurement_length = 31556952 / 12
t = np.arange(0.0, measurement_length, sample_step_size)
target_signal = np.sin(f_target * 2 * np.pi * (t - delta_target)) * np.exp(-b_target * (t - delta_target) ** 2 / 2) * a_target
# signal=1.01**t*np.sin(0.15*t**(1/2)*np.pi*(t-70))*np.exp(-0.5*(t-70)**2/25)*30
template = lambda f, b, a, delta, t: np.sin(f * 2 * np.pi * (t - delta)) * np.exp(-b * (t - delta) ** 2 / 2) * a
# template = lambda a, b, c, delta, t: a**t*np.sin(b*t**(1/2)*np.pi*(t-delta))*np.exp(-0.5*(t-delta)**2/25)*c
# def template(a, f, c, delta, t):
#     if t < delta:
#         amplitude = a**t
#     else:
#         amplitude = a**(delta-t)
#     return amplitude*np.sin(f*t*np.pi*(t-delta))*c


# Make some noise. ------------------------------------------------
# Plain random noise.
# Set the noise level.
SNR = 1
noise_level = a_target / SNR
n1 = np.zeros(len(target_signal))
n1[1:-1] = np.random.normal(0, 1, len(target_signal) - 2)

# Smoothing.
for it in range(1):
    n1_new = n1
    n1_new[1:-1] = (n1[0:-2] + n1[2:] + n1[1:-1]) / 3.0
# n1_new = n1
n1 = n1_new / np.max(np.abs(n1_new))

# Add scaled version to the time series.
data = target_signal + n1 * noise_level

sig_fft = sp.fft.fft(data)
power = np.abs(sig_fft)
sample_freq = sp.fft.fftfreq(data.size, d=sample_step_size)
pos_mask = np.where(sample_freq > 0)
freqs = sample_freq[pos_mask]
peak_freq = freqs[power[pos_mask].argmax()]
high_freq_fft = sig_fft.copy()
high_freq_fft[np.abs(sample_freq) > peak_freq * 1.6] = 0
filtered_fft = high_freq_fft.copy()
filtered_fft[np.abs(sample_freq) < peak_freq * 0.4] = 0
data_filtered = sp.fft.ifft(filtered_fft)

f, Pxx_den = sp.signal.periodogram(data, sample_rate)
plt.figure()
plt.semilogy(f, Pxx_den, label='power noisy')
plt.semilogy(sample_freq, power, label='power')
plt.axvline(f_target, color='r')
plt.xlabel('frequency [Hz]')
plt.ylabel('PSD [V**2/Hz]')
plt.legend()
plt.show()

plt.figure()
plt.plot(t, data, label='Noisy signal', color='k', linewidth=0.5)
plt.plot(t, target_signal, 'r', label='True signal')
plt.plot(t, data_filtered, label='Filtered signal', color='b')
plt.xlabel('time (seconds)')
plt.grid(True)
plt.axis('tight')
plt.legend(loc='upper right')

# Sample rate and desired cutoff frequencies (in Hz).
lowcut = 0.02
highcut = 0.07


# Apply butter bandpass filter
# data_filtered = butter_bandpass_filter(data, lowcut, highcut, sample_rate, order=10)
# plt.figure()
# plt.plot(t, data, label='Noisy signal', color = 'k', linewidth=0.5)
# plt.plot(t,target_signal,'r', label='True signal')
# plt.plot(t, data_filtered, label='Filtered signal', color='b')
# plt.xlabel('time (seconds)')
# plt.grid(True)
# plt.axis('tight')
# plt.legend(loc='upper right')
# Define the prior in data space. ----------------------------------

def prior_data(data, s2, noise_level):
    p = np.sum((data - s2) ** 2) / (len(data) * noise_level ** 2)
    return np.exp(-p / 2.0)


# Some input. ------------------------------------------------------

# Number of histogram bins.
n_bin = 50
# Total number of proposed samples.
number_of_samples = 10 ** 3
number_of_parameters = 4

medium_jump = 5
large_jump = 10

# Make the first random sample. ------------------------------------

samples = np.zeros((number_of_samples, number_of_parameters))

# Random time shift between 10 and 30.
samples[0, 0] = f_target * 10 * np.random.rand()
samples[0, 1] = 0.01 * np.random.rand()
samples[0, 2] = 100.0 * np.random.rand()
samples[0, 3] = measurement_length * np.random.rand()

# Make test time series by shifting s2.
s = template(samples[0, 0], samples[0, 1], samples[0, 2], samples[0, 3], t)

# Evaluate posterior for the first sample.
p = prior_data(s, data, noise_level)

# Run Metropolis-Hastings sampler. ---------------------------------

for i in range(1, number_of_samples):

    # Random time shift proposal.
    # large jump
    if np.mod(i, large_jump) == 0:
        f_test = f_target * 10 * np.random.rand()
        b_test = 0.01 * np.random.rand()
        a_test = 100 * np.random.rand()
        delta_test = measurement_length * np.random.rand()
    # medium jump
    elif np.mod(i, medium_jump) == 0:
        f_test = samples[i - 1, 0] + f_target * (np.random.rand() - 0.5)
        b_test = samples[i - 1, 1] + 0.001 * (np.random.rand() - 0.5)
        a_test = samples[i - 1, 2] + 10 * (np.random.rand() - 0.5)
        delta_test = samples[i - 1, 3] + measurement_length / 10 * (np.random.rand() - 0.5)
    # small jump
    else:
        f_test = samples[i - 1, 0] + f_target / 10 * (np.random.rand() - 0.5)
        b_test = samples[i - 1, 1] + 0.0001 * (np.random.rand() - 0.5)
        a_test = samples[i - 1, 2] + 1 * (np.random.rand() - 0.5)
        delta_test = samples[i - 1, 3] + measurement_length / 10 * (np.random.rand() - 0.5)

    if i > number_of_samples / 5:
        f_test = samples[i - 1, 0]
    # a_test = 0.5
    # b_test = 0
    # delta_test = 70

    # Make test time series by shifting s2.
    s = template(f_test, b_test, a_test, delta_test, t)

    # Evaluate posterior.
    p_test = prior_data(s, data_filtered, noise_level)
    T_inv = i / number_of_samples * 10 ** 7
    # Apply Metropolis rule.
    if (p_test / p) ** T_inv > np.random.rand():  # L^i/L^j
        p = p_test
        samples[i, 0] = f_test
        samples[i, 1] = b_test
        samples[i, 2] = a_test
        samples[i, 3] = delta_test
    else:
        samples[i] = samples[i - 1]

# Plot results. ----------------------------------------------------
# %%
n_bin = 200
plt.figure()
plt.suptitle("sampled posterior")
plt.subplot(221)
plt.axvline(x=f_target, color='r')
n, bins, patches = plt.hist(samples[:, 0], n_bin, density=True, facecolor='k', alpha=0.7)
plt.xlabel('frequency')

plt.subplot(222)
plt.axvline(x=b_target, color='r')
n, bins, patches = plt.hist(samples[:, 1], n_bin, density=True, facecolor='k', alpha=0.5)
plt.xlabel('b')

plt.subplot(223)
plt.axvline(x=a_target, color='r')
n, bins, patches = plt.hist(samples[:, 2], n_bin, density=True, facecolor='k', alpha=0.5)
plt.xlabel('amplitude')

plt.subplot(224)
plt.axvline(x=delta_target, color='r')
n, bins, patches = plt.hist(samples[:, 3], n_bin, density=True, facecolor='k', alpha=0.5)
plt.xlabel('time shift [s]')
plt.savefig(r'pictures\time_histo.png')
plt.tight_layout

plt.figure()
plt.subplot(221)
plt.axvline(x=f_target, color='r')
plt.plot(samples[:, 0], range(number_of_samples), 'k')
plt.xlabel('frequency [Hz]')
plt.ylabel('sample number')
plt.subplot(222)
plt.axvline(x=b_target, color='r')
plt.plot(samples[:, 1], range(number_of_samples), 'k')
plt.xlabel('b')
plt.subplot(223)
plt.axvline(x=a_target, color='r')
plt.plot(samples[:, 2], range(number_of_samples), 'k')
plt.xlabel('amplitude')
plt.ylabel('sample number')
plt.subplot(224)
plt.axvline(x=delta_target, color='r')
plt.plot(samples[:, 3], range(number_of_samples), 'k')
plt.xlabel('time shift [s]')
plt.savefig(r'pictures\time_samples.png')

plt.figure()
s = template(samples[-1, 0], samples[-1, 1], samples[-1, 2], samples[-1, 3], t)
plt.plot(t, data, label='Noisy signal', color='k')
plt.plot(t, data_filtered, 'b', label='Filtered signal')
plt.plot(t, target_signal, 'r', linewidth=1, label='True signal')
plt.plot(t, s, 'g', linewidth=1, alpha=0.9, label='Found signal')
plt.xlabel('t [s]')
plt.legend()
plt.title('time series')
plt.savefig(r'pictures\time series.png')
plt.show()

print(samples[-1, :])

# %%
