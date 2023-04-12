""" Provide time domain window functions.
"""

import numpy as np
from scipy import signal

#pylint:disable=C0103

def window(tm, margin=1000.0, kap=0.005, show=False, xl=None):
    """Return time domain window function to taper a margin in
    seconds at both ends.

    Args:
        tm (array_like): time axis
        margin (float): margin length to taper
        kap (float): control how abrupt the window edges are
        show (bool): if True, the window will be plotted
        xl (float): margin length to taper. *Deprecated*

    """
    if xl is not None:
        margin = xl
        import warnings
        warnings.warn("xl is deprecated; use margin", DeprecationWarning)
    xr = tm[-1] - margin
    xl = tm[0] + margin
    winl = 0.5*(1.0 + np.tanh(kap*(tm-xl)))
    winr = 0.5*(1.0 - np.tanh(kap*(tm-xr)))
    if show:
        import matplotlib.pyplot as plt
        plt.plot(tm, winl)
        plt.plot(tm, winr)
        plt.grid(True)
        plt.show()
    return (winl*winr)

def tukey(N, alpha):
    """ Return time domain window function of size N.
    """
    # alpha -- parameter the defines the shape of the window
    w         = np.zeros(N, dtype=float)
    i         = np.arange(0,N,1)
    r         = (2.0*i)/(alpha*(N-1))
    l1        = int(np.rint((alpha*(N-1))/2.0))
    l2        = int(np.rint((N-1)*(1.0-alpha/2.0)))
    w[0:l1]   = 0.5*(1.0 + np.cos(np.pi*(r[0:l1] - 1.0)))
    w[l1:l2]  = 1.0
    w[l2:N-1] = 0.5*(1+np.cos(np.pi*(r[l2:N-1] - 2.0/alpha + 1.0)))
    return w

def butter_lowpass_filter(data, cutoff, fs, order=4, show=False):
    """ Return low pass filtered data
    """
    #nyq = 0.5 * fs  # Nyquist Frequency
    #normal_cutoff = cutoff / nyq
    # Get the filter coefficients
    b, a = signal.butter(order, cutoff, btype='low', analog=False)
    y = signal.filtfilt(b, a, data)

    if show:
        import matplotlib.pyplot as plt
        w, h = signal.freqz(b, a)
        plt.figure()
        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    return y


def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):
    """ Return low pass filtered data. 
    """
    nyq = 0.5 * fs
    low, high = lowcut / nyq, highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    y = signal.filtfilt(b, a, data)
    return y

def window_planck_vec(x, xi, xf, deltaxi, deltaxf):
    """Definitions for the windowing function (from lisabeta). 
      
    In order to avoid overflows in the exponentials, we set the
    boundaries (di, df) so that anything below 10^-20 is considered
    zero
    """

    di = deltaxi/(20*np.log(10))
    df = deltaxf/(20*np.log(10))
    w = np.zeros(len(x), dtype=float)
    mask = np.logical_or(x <= xi + di, x >= xf - df)
    w[mask] = 0.
    mask = np.logical_and(xi + di < x, x < xi + deltaxi - di)
    xm = x[mask]
    w[mask] = 1./(1 + np.exp(deltaxi/(xm - xi) + deltaxi/(xm - (xi + deltaxi))))
    mask = np.logical_and(xi + deltaxi - di <= x, x <= xf - deltaxf + df)
    w[mask] = 1.
    mask = np.logical_and(xf - deltaxf + df < x, x < xf - df)
    xm = x[mask]
    w[mask] = 1./(1 + np.exp(-(deltaxf/(xm - (xf - deltaxf))) - deltaxf/(xm - xf)))
    return w
