import gwpy
from gwosc.datasets import event_gps
from gwpy.timeseries import TimeSeries
import pylab
import pycbc.psd
from pycbc import frame
from pycbc.filter import resample_to_delta_t, highpass, matched_filter
from pycbc.waveform import get_td_waveform
import matplotlib as plt

# An example of how to read the data from these files:
file_name = "challenge3.gwf"

# LOSC bulk data typically uses the same convention for internal channels names
# Strain is typically IFO:LOSC-STRAIN, where IFO can be H1/L1/V1.
channel_name = "H1:CHALLENGE3"

start = 0
end = start + 128

hdata = frame.read_frame(file_name, channel_name)

print("GW170817 data")
print(hdata)
# Remove the low frequency content and downsample the data to 2048Hz
hdata = highpass(hdata, 15.0)
hdata = resample_to_delta_t(hdata, 1.0/2048)

pylab.plot(hdata.sample_times, hdata)
pylab.xlabel('Time (s)')
pylab.show()

sample_rate = hdata.sample_rate
data_length = hdata.duration
delta_t = 1.0 / sample_rate
flow = 10
delta_f = 1.0 / data_length
flen = int(sample_rate / (2 * delta_f)) + 1
psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)

seg_len = int(4 / delta_t)
seg_stride = int(seg_len / 2)
estimated_psd = pycbc.psd.welch(hdata,seg_len=seg_len,seg_stride=seg_stride)

pylab.loglog(estimated_psd.sample_frequencies, estimated_psd, label='data')
pylab.loglog(psd.sample_frequencies, psd, linewidth=3, label='known psd')
pylab.xlim(xmin=20, xmax=sample_rate/2)
pylab.ylim(1e-47, 1e-45)
pylab.xlabel('Frequency [Hz]')
pylab.ylabel('Power spectral density')
pylab.legend()
pylab.grid()
pylab.show()

# times, freqs, power = hdata.qtransform(.001, logfsteps=100, qrange=(8, 8), frange=(20, 512))
# pylab.figure(figsize=[15, 3])
# pylab.pcolormesh(times, freqs, power**0.5)
# # pylab.xlim(xmin=-15, xmax=-14.5)
# pylab.yscale('log')
# pylab.show()

conditioned = hdata.crop(2, 2)
hp, _ = get_td_waveform(approximant='SEOBNRv4_opt',
                         mass1=10,
                         mass2=10,
                         delta_t=1.0/sample_rate,
                         f_lower=25)
hp.resize(len(conditioned))
template = hp.cyclic_time_shift(hp.start_time)

# We use 4 second samples of our time series in Welch method.
psd = conditioned.psd(4)

# Now that we have the psd we need to interpolate it to match our data
# and then limit the filter length of 1 / PSD. After this, we can
# directly use this PSD to filter the data in a controlled manner
psd = pycbc.psd.interpolate(psd, conditioned.delta_f)

# 1/PSD will now act as a filter with an effective length of 4 seconds
# Since the data has been highpassed above 15 Hz, and will have low values
# below this we need to inform the function to not include frequencies
# below this frequency. 
psd = pycbc.psd.inverse_spectrum_truncation(psd, int(4 * conditioned.sample_rate),
                                  low_frequency_cutoff=15)

snr = matched_filter(template, conditioned,
                     psd=psd, low_frequency_cutoff=20)
snr = snr.crop(4 + 4, 4)

# Why are we taking an abs() here?
# The `matched_filter` function actually returns a 'complex' SNR.
# What that means is that the real portion correponds to the SNR
# associated with directly filtering the template with the data.
# The imaginary portion corresponds to filtering with a template that
# is 90 degrees out of phase. Since the phase of a signal may be 
# anything, we choose to maximize over the phase of the signal.
pylab.figure(figsize=[10, 4])
pylab.plot(snr.sample_times, abs(snr))
pylab.ylabel('Signal-to-noise')
pylab.xlabel('Time (s)')
pylab.show()

peak = abs(snr).numpy().argmax()
snrp = snr[peak]
time = snr.sample_times[peak]

print("We found a signal at {}s with SNR {}".format(time, 
                                                    abs(snrp)))