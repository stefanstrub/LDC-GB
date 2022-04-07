import gwpy
from gwosc.datasets import event_gps
from gwpy.timeseries import TimeSeries
import pylab
import pycbc.psd
from pycbc import frame
from pycbc.filter import resample_to_delta_t, highpass
import matplotlib as plt

# An example of how to read the data from these files:
DATAPATH = "/home/stefan/GW_open/"
file_name = "challenge1.gwf"

# LOSC bulk data typically uses the same convention for internal channels names
# Strain is typically IFO:LOSC-STRAIN, where IFO can be H1/L1/V1.
channel_name = "H1:CHALLENGE1"

start = 0
end = start + 128

hdata = frame.read_frame(file_name, 'H1:CHALLENGE1')

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
delta_f = 1.0 / 128
flen = int(sample_rate / (2 * delta_f)) + 1
psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)

seg_len = int(4 / delta_t)
seg_stride = int(seg_len / 2)
estimated_psd = pycbc.psd.welch(hdata,seg_len=seg_len,seg_stride=seg_stride)

pylab.loglog(estimated_psd.sample_frequencies, estimated_psd, label='data')
pylab.loglog(psd.sample_frequencies, psd, linewidth=3, label='known psd')
pylab.xlim(xmin=flow, xmax=512)
pylab.ylim(1e-47, 1e-45)
pylab.ylim(1e-47, 1e-45)
pylab.xlabel('Frequency [Hz]')
pylab.ylabel('Power spectral density')
pylab.legend()
pylab.grid()
pylab.show()

times, freqs, power = hdata.qtransform(.001, logfsteps=100, qrange=(8, 8), frange=(20, 512))
pylab.figure(figsize=[15, 3])
pylab.pcolormesh(times, freqs, power**0.5)
pylab.yscale('log')
pylab.show()

print('end')