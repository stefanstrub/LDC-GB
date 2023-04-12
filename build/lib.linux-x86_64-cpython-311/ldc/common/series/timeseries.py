"""
Time and frequency series containers based on xarray.
"""
import datetime
import xarray as xr
import numpy as np
import astropy.units as u

#pylint:disable=C0103

def TimeSeries(array, t0=None, dt=None, units=None, name=None, ts=None):
    """Represent an equispaced LDC time series by way of an xarray.

    Args:
        array (numpy-like object): the data
        t0 (double or numpy.datetime64): time of first sample (in seconds if given as double)
        dt (double or numpy.timedelta64): sample cadence (in seconds if given as double)
        units (str): physical units (MKS + LDC additions)
        ts (array_like): time axis. If provided, `dt` and `t0` will be ignored.

    Returns:
        xarray.DataArray: the wrapped and annotated data

    TODO:
        accept more meta data

    >>> t = np.linspace(0,8,1000)
    >>> s = np.sin(2 * np.pi * (t + 0.1 * t**2))
    >>> hp = TimeSeries(s, t0=datetime.datetime.fromisoformat('2030-01-01 00:00:00'), dt=10, name='hp', units='strain')
    """



    if t0 is not None and isinstance(t0, datetime.datetime):
        t0 = t0.timestamp()

    # build the time axis
    if ts is None:
        ts = xr.DataArray(np.arange(t0, t0 + len(array) * dt, dt), dims=('t'))
    else:
        assert t0 is not None and dt is not None, "In absence of `ts`, one should provide both `t0` and `dt`."
        ts = xr.DataArray(ts, dims=('t'))

    # if the time axis is numeric, give it units of second
    if not np.issubdtype(ts.dtype, np.datetime64):
        ts.attrs['units'] = u.s

    return xr.DataArray(array, dims=('t'), coords={'t': ts},
                        name=name, attrs={'units': units, 't0': t0, 'dt': dt})


def FrequencySeries(array, df=None, kmin=None, t0=0, units=None, name=None, fs=None):
    """Represent a frequency-domain LDC signal by way of an xarray.

    Args:
        array (numpy-like object): the data
        df (double): frequency spacing (in Hertz)
        kmin: index of the first sample (fmin = kmin * df)
        t0 (double or numpy.datetime64): time offset of the signal (in seconds if given as double)
        units (str): physical units of the corresponding TimeSeries (MKS + LDC additions)
        fs (array_like): frequency axis. If provided, `df` and `kmin` will be ignored.

    Returns:
        xarray.DataArray: the wrapped and annotated data

    TODO:
        accept more metadata
        use the actual FFT units

    >>> t = np.linspace(0,8,1000)
    >>> s = np.sin(2 * np.pi * (t + 0.1 * t**2))
    >>> hp = TimeSeries(s, name='hp', units='strain', t0=0, dt=1)
    >>> hp.ts.fft().attrs
    {'units': 'strain', 'df': 0.001, 'kmin': 0, 't0': 0}
    """
    # build the frequency axis or check the one that was provided
    if fs is None:
        assert kmin is not None and df is not None, "In absence of `fs`, one should provide both `kmin` and `df`."
        fmin = kmin * df
        fs = xr.DataArray(df * np.arange(kmin, kmin + len(array)), dims=('f'))
    else:
        fmin = fs[0]
        df = fs[1] - fs[0]
        kmin = int(np.rint(fmin/df))
        assert np.isclose(kmin*df, fmin), f"The `fs` provided is not correct: {kmin}*{df}=={kmin*df}!={fmin}."
        fs = xr.DataArray(fs, dims=('f'))
    # add units to frequency axis
    fs.attrs['units'] = u.Hz
    return xr.DataArray(array, dims=('f'), coords={'f': fs},
                    name=name, attrs={'units': units, 'df': df, 'kmin': kmin, 't0': t0})

@xr.register_dataarray_accessor('ts')
class TimeSeriesAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    def fft(self, win=None, **kwargs):
        """Obtain the frequency-domain FrequencySeries representation of a
        TimeSeries.

        >>> t = np.linspace(0,8,1000)
        >>> s = np.sin(2 * np.pi * (t + 0.1 * t**2))
        >>> hp = TimeSeries(s, t0=0, dt=1, name='hp', units='strain')
        >>> hp.ts.fft().f.values[0:5]
        array([0.   , 0.001, 0.002, 0.003, 0.004])
        """
        array = self._obj
        name = array.name
        t0, dt, units = array.attrs['t0'], array.attrs['dt'], array.attrs['units']
        times = array['t'].values
        if win is not None:
            array = array.copy()*win(times, **kwargs)
        frequencies = np.fft.rfftfreq(array.size, d=dt)
        return FrequencySeries(np.fft.rfft(array)*dt,
                               fs=frequencies, t0=t0,
                               units=units, name=name)

    def ifft(self, dt=None, win=None, **kwargs):
        """Obtain the real-space TimeSeries representation of a FrequencySeries.

        Args:
        dt (double): TimeSeries cadence (defaults to 0.5 / highest f represented in FrequencySeries)

        >>> t = np.linspace(0,8,1000)
        >>> s = np.sin(2 * np.pi * (t + 0.1 * t**2))
        >>> hp = TimeSeries(s, t0=0, dt=10)
        >>> hpf = hp.ts.fft()
        >>> np.allclose(hpf.ts.ifft().values, hp.values)
        True
        """

        array = self._obj
        name = array.name
        df = array.attrs['df']
        t0, units= array.attrs['t0'], array.attrs['units']
        #Finding the default time step assumed there
        frequencies = array['f'].values
        highest_f = frequencies[-1]
        if dt is None:
            dt = 0.5 / highest_f
        # Allow tapering before irfft
        if win is not None:
            array = array.copy()*win(frequencies, **kwargs)

        target_length = int(1 / df / dt)

        kmin = array.attrs['kmin']
        # The first element of array should correspond to frequency zero
        array = np.pad(array, (kmin, 0))
        return TimeSeries(np.fft.irfft(array / dt, n=target_length),
                          t0=t0, dt=dt,
                          units=units, name=name)

if __name__ == "__main__":

    import doctest
    doctest.testmod()
