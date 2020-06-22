""" High level waveforms and TDI functions. 
"""
import numpy as np
import ldc.waveform.fastGB as FB
from .hphc import HpHc

def get_td_waveform(delta_t, start_time, duration, source_type,
                    approximant, name="", **kwargs):
    """ Return hp,hc in time domain. 
    """
    
    GW = HpHc.type(name, source_type, approximant)
    GW.set_param(kwargs)
    time_vector = np.arange(start_time, duration, delta_t)
    return GW.compute_hphc_td(time_vector)
    
    

def get_fd_tdixyz(delta_t, duration, source_type, approximant, **kwargs):
    """Return tdi x,y,z in Fourier domain

    TODO: as for h+hc, use an intermediate abstract class to avoid
    duplication of if/raise statements
    """

    if source_type=="GB":
        GB = FB.FastGB(delta_t=delta_t, T=duration) # in seconds
        freqT, X, Y, Z = GB.get_fd_tdixyz(template=kwargs)
        return freqT, X, Y, Z
    else:
        raise NotImplementedError
    
