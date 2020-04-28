
# HpHx

## Proposal

Proposal: implement the pycbc.waveform.waveform functional interface,
i.e.

```
hp, hc = get_td_waveform(template=None, **kwargs)
```

Here kwargs includes all the source parameters (in fixed units and
with standardized names) and the approximant (a string). The sampling
cadence is given as `delta_t` and the length is determined by the
initial frequency `f_lower`, as appropriate for a coalescence, but in
a LISA context we may instead specify the initial epoch and
duration. These parameters can also be taken from template, ducktyped
as a namespace. The returned hp and hc are
pycbc.types.timeseries.TimeSeries objects, which are ndarray-like (or
PyCUDA-like!) but include extra information such as `delta_t`,
`start_time`, and methods implementing Fourier transforms, and
more. Waveforms are (apparently) aligned at the nominal coalescence
time.

```
fp, fc = get_fd_waveform(template=None, **kwargs)
```

Same, with sampling parameters `delta_f` and `f_lower`. Outputs
pycbc.types.frequencyseries.FrequencySeries objects, again with smart
methods.

## Current status

Only time domain is available for now. See waveform/test/demo.py for
an example. The ouput hp,hc are simple numpy array, not encapsulated
into pycbc timeseries.  


# Fast approximation

## Proposal

For classes that implement fast TDI approximants, we will also define

```
x, y, z = get_td_tdixyz(template=None, **kwargs)
fx, fy, fz = get_fd_tdixyz(template=None, **kwargs)
```

Here kwargs includes also orbits (the LDC object) and tdigeneration
(and perhaps tdiobservable). The LDC TDI approximants that hardwire
the orbits should check that they match the requested orbit. The
output can be a derived class TDITimeSeries that includes orbit and
TDI information.

## Current status

Only fastGB available for now. See waveform/test/demo.py for an
example. An abstract object like for hphx is still missing. The ouput
x,y,z are simple numpy array, not encapsulated into pycbc
freq. series.



