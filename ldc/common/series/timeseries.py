""" 
Time and frequency series containers based on xarray. 
"""
import xarray as xr
import numpy as np


def TimeSeries(array, t0=0, dt=1, units=None, name=None, ts=None):
    """Represent an equispaced LDC time series by way of an xarray.

    Args:
        array (numpy-like object): the data
        t0 (double or numpy.datetime64): time of first sample (in seconds if given as double)
        dt (double or numpy.timedelta64): sample cadence (in seconds if given as double)
        units (str): physical units (MKS + LDC additions) 
        
    Returns:
        xarray.DataArray: the wrapped and annotated data
        
    TODO:
        accept more meta data
    """

    # build the time axis or check the one that was provided
    if ts is None:
        ts = xr.DataArray(np.arange(t0, t0 + len(array) * dt, dt), dims=('t'))
    else:
        assert (ts[0] == t0) and (ts[1] - ts[0] == dt)
        ts = xr.DataArray(ts, dims=('t'))
    
    # if the time axis is numeric, give it units of second
    if not np.issubdtype(ts.dtype, np.datetime64):
        ts.attrs['units'] = 's'
        
    return xr.DataArray(array, dims=('t'), coords={'t': ts},
                        name=name, attrs={'units': units, 't0': t0, 'dt': dt})


def FrequencySeries(array, df=0, kmin=0, t0=0, units=None, name=None):
    """Represent a frequency-domain LDC signal by way of an xarray.

    Args:
        array (numpy-like object): the data
        df (double): frequency spacing (in Hertz)
        kmin: the index of the first sample (so that its frequency is (kmin * df) Hertz)
        t0 (double or numpy.datetime64): time offset of the signal (in seconds if given as double)
        units (str): physical units of the corresponding TimeSeries (MKS + LDC additions)
        
    Returns:
        xarray.DataArray: the wrapped and annotated data
        
    TODO:
        accept more metadata
        use the actual FFT units
    """

    # build the frequency axis
    fs = xr.DataArray(df * np.arange(kmin, kmin + len(array)), dims=('f'), attrs={'units': '1/s'})
            
    return xr.DataArray(array, dims=('f'), coords={'f': fs},
                        name=name, attrs={'units': units, 'df': df, 'kmin': kmin, 't0': t0})

@xr.register_dataarray_accessor('ts')
class TimeSeriesAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
 
    def fft(self):
        """Obtain the frequency-domain FrequencySeries representation of a TimeSeries.        
        """
        
        array = self._obj
        t0, dt, units = array.attrs['t0'], array.attrs['dt'], array.attrs['units']     
        
        return FrequencySeries(np.fft.rfft(array)*dt,
                               df=1.0/(dt * len(array)), kmin=0, t0=t0,
                               units=units, name=array.name)
    
    def ifft(self, dt=None):
        """Obtain the real-space TimeSeries representation of a FrequencySeries.
                
        Args:
            dt (double): TimeSeries cadence (defaults to 0.5 / highest f represented in FrequencySeries)
        """
        
        array = self._obj
        df, kmin, t0, units = array.attrs['df'], array.attrs['kmin'], array.attrs['t0'], array.attrs['units']

        if dt is None:
            n = 2 * (kmin + len(array) - 1)
            dt = 1 / df / n
        else:
            #n = int(1.0 / (dt * df))
            n = int(1.0 / (dt * df))

            
        # pad the array if needed
        padded = np.zeros(int(n/2 + 1), dtype=array.dtype)
        padded[kmin:(kmin + len(array))] = array.values[:]
        padded *= df*n
        
        return TimeSeries(np.fft.irfft(padded),
                          t0=t0, dt=dt,
                          units=units, name=array.name)


if __name__ == "__main__":
    
    import doctest
    t = np.linspace(0,8,1000)
    s = np.sin(2 * np.pi * (t + 0.1 * t**2))
    c = np.cos(2 * np.pi * (t + 0.1 * t**2))
    hp = TimeSeries(s, name='hp', units='strain')
    hp_f = hp.ts.fft()
    #hpf = fhp.ts.ifft()
    
    
