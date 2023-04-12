""" High level waveforms and TDI functions.
"""
import numpy as np

# the second form does not work in Python 3.6
try:
    from ldc.waveform import fastGB as FB
except ImportError:
    FB = None

# import ldc.waveform.fastGB as FB

from ldc.common.series import TimeSeries
from .hphc import HpHc

#pylint:disable=C0103

def get_td_waveform(delta_t, start_time, duration, source_type,
                    approximant, name="", **kwargs):
    """ Return hp,hc in time domain as TimeSeries object.
    """

    GW = HpHc.type(name, source_type, approximant)
    GW.set_param(kwargs)
    time_vector = np.arange(start_time, duration, delta_t)
    hp, hc = GW.compute_hphc_td(time_vector)
    return (TimeSeries(hp, name="hp", units='strain', dt=delta_t, t0=start_time),
            TimeSeries(hc, name="hc", units='strain', dt=delta_t, t0=start_time))



def get_fd_tdixyz(delta_t, duration, source_type, approximant, **kwargs):
    """Return tdi x,y,z in Fourier domain as FrequencySeries object.

    TODO: as for h+hc, use an intermediate abstract class to avoid
    duplication of if/raise statements
    """

    if source_type=="GB" and FB is not None:
        GB = FB.FastGB(delta_t=delta_t, T=duration) # in seconds
        X, Y, Z = GB.get_fd_tdixyz(template=kwargs)
        return X, Y, Z
    else:
        if source_type=="GB":
            print("FastGB not found, please install without --no-fastGB")
        raise NotImplementedError
