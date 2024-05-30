""" Provide time domain window functions.
"""

import numpy as np
from scipy import signal
from scipy.signal import windows
from scipy.special import expit
from functools import wraps
#pylint:disable=C0103

# A class might be more suitable
def _add_show(name):
    """Decorator adding illustration codes to windows."""
    def decorator(func):
        @wraps(func)
        def decorated_func(time_axis, *args, show=False, **kwargs):
            _window = func(time_axis, *args, **kwargs)
            if show:
                import matplotlib.pyplot as plt
                fig, axs = plt.subplots(3)
                fig.tight_layout(h_pad=4)
                # Plot window in TD
                axs[0].plot(time_axis, _window)
                axs[0].set_title(f"{name} window")
                axs[0].set_xlabel("Time grid")
                axs[0].set_ylabel("Amplitude")
                axs[0].set_ylim([0, 1.1])
                # We want enough resolution in frequency
                if _window.shape[0] > 2048:
                    N = _window.shape[0]
                else:
                    N = 2048
                _fd_window = np.fft.fft(_window, N)
                freq = np.fft.fftfreq(N, d=time_axis[1]-time_axis[0])
                response = 20 * np.log10(np.abs(_fd_window / np.abs(_fd_window).max()))
                # Plot frequency response
                axs[1].plot(np.fft.fftshift(freq), np.fft.fftshift(response))
                axs[1].set_title(f"Frequency response of the {name} window")
                axs[1].set_ylabel("Normalized magnitude [dB]")
                axs[1].set_xlabel("Frequency [Hz]")
                axs[1].set_ylim([-120, 0])
                # Same but zoomed in frequency
                axs[2].plot(np.fft.fftshift(freq), np.fft.fftshift(response))
                axs[2].set_title(f"Frequency response of the {name} window (zoomed in)")
                axs[2].set_ylabel("Normalized magnitude [dB]")
                axs[2].set_xlabel("Frequency [Hz]")
                axs[2].set_ylim([-120, 0])
                axs[2].set_xlim([-0.001, 0.001])
                return fig, axs, _window
            return _window
        return decorated_func
    return decorator
    

@_add_show(name="Default")
def window(tm, margin=1000.0, kap=0.005, xl=None):
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
    return (winl*winr)

@_add_show(name="Tukey")
def tukey(time_axis, alpha=0.5, sym=False):
    """ Return time domain Tukey window function, also known as tapered cosine window.
    
    Args:
        time_axis (array_like): Time axis
        alpha (float): Shape parameter of the Tukey window, representing the fraction of the window inside the cosine tapered region. 
            If zero, the Tukey window is equivalent to a rectangular window. 
            If one, the Tukey window is equivalent to a Hann window.
        sym (bool): When False (default), generates a periodic window, for use in spectral analysis.
            When True, generates a symmetric window, for use in filter design. 
        show (bool): If True, the window and its frequency response will be plotted
        
    """
    _window = windows.tukey(time_axis.size, alpha=alpha, sym=sym)
    return _window

@_add_show(name="flattop")
def flattop(time_axis, sym=False):
    """Return a flat top window.

    Parameters
    ----------
    time_axis (array_like): Time axis
    sym : bool, optional
        When True, generates a symmetric window, for use in filter
        design.
        When False (default), generates a periodic window, for use in spectral analysis.

    Returns
    -------
    w : ndarray
        The window, with the maximum value normalized to 1 (though the value 1
        does not appear if the length of `time_axis` is even and `sym` is True).

    """
    _window = windows.flattop(M=time_axis.size, sym=sym)
    return _window

@_add_show(name="General cosine")
def general_cosine(time_axis, a, sym=False):
    r"""
    Generic weighted sum of cosine terms window

    Parameters
    ----------
    time_axis (array_like): Time axis
    a : array_like
        Sequence of weighting coefficients. This uses the convention of being
        centered on the origin, so these will typically all be positive
        numbers, not alternating sign.
    sym : bool, optional
        When True (default), generates a symmetric window, for use in filter
        design.
        When False, generates a periodic window, for use in spectral analysis.

    Returns
    -------
    w : ndarray
        The array of window values.

    References
    ----------
    .. [1] A. Nuttall, "Some windows with very good sidelobe behavior," IEEE
           Transactions on Acoustics, Speech, and Signal Processing, vol. 29,
           no. 1, pp. 84-91, Feb 1981. :doi:`10.1109/TASSP.1981.1163506`.
    .. [2] Heinzel G. et al., "Spectrum and spectral density estimation by the
           Discrete Fourier transform (DFT), including a comprehensive list of
           window functions and some new flat-top windows", February 15, 2002
           https://holometer.fnal.gov/GH_FFT.pdf

    Examples
    --------
    Heinzel describes a flat-top window named "HFT90D" with formula: [2]_

    .. math::  w_j = 1 - 1.942604 \cos(z) + 1.340318 \cos(2z)
               - 0.440811 \cos(3z) + 0.043097 \cos(4z)

    where

    .. math::  z = \frac{2 \pi j}{N}, j = 0...N - 1

    Since this uses the convention of starting at the origin, to reproduce the
    window, we need to convert every other coefficient to a positive number:

    >>> HFT90D = [1, 1.942604, 1.340318, 0.440811, 0.043097]
    """
    _window = windows.general_cosine(M=time_axis.size, a=a, sym=sym)
    return _window

@_add_show(name="Planck")
def planck(time_axis, left_margin=0, right_margin=0):
    """ Return time domain Planck window function. 
    Adapted from gwpy (https://gwpy.github.io/).
    
    Args:
        time_axis (array_like): Time axis
        left_margin (float): Margin to taper on the left, should be less than half width
        right_margin (float): Margin to taper on the right, should be less than half width
        show (bool): If True, the window and its frequency response will be plotted 
    
    References
    ----------
    .. [1] McKechan, D.J.A., Robinson, C., and Sathyaprakash, B.S. (April
        2010). "A tapering window for time-domain templates and simulated
        signals in the detection of gravitational waves from coalescing
        compact binaries". Classical and Quantum Gravity 27 (8).
        :doi:`10.1088/0264-9381/27/8/084020`
    .. [2] Wikipedia, "Window function",
        https://en.wikipedia.org/wiki/Window_function#Planck-taper_window
    """
    N = time_axis.shape[0]
    cadence = time_axis[1] - time_axis[0]
    nleft = int(np.rint(left_margin / cadence))
    nright = int(np.rint(right_margin / cadence))

    # This is copied from gwpy
    # -- Planck taper window ------------------------------------------------------
    # source: https://arxiv.org/abs/1003.2939
    def gwpy_planck(N, nleft=0, nright=0):
        """Return a Planck taper window.
        Parameters
        ----------
        N : `int`
            Number of samples in the output window
        nleft : `int`, optional
            Number of samples to taper on the left, should be less than `N/2`
        nright : `int`, optional
            Number of samples to taper on the right, should be less than `N/2`
        Returns
        -------
        w : `ndarray`
            The window, with the maximum value normalized to 1 and at least one
            end tapered smoothly to 0.
        Examples
        --------
        To taper 0.1 seconds on both ends of one second of data sampled at 2048 Hz:
        >>> from gwpy.signal.window import planck
        >>> w = planck(2048, nleft=205, nright=205)
        References
        ----------
        .. [1] McKechan, D.J.A., Robinson, C., and Sathyaprakash, B.S. (April
            2010). "A tapering window for time-domain templates and simulated
            signals in the detection of gravitational waves from coalescing
            compact binaries". Classical and Quantum Gravity 27 (8).
            :doi:`10.1088/0264-9381/27/8/084020`
        .. [2] Wikipedia, "Window function",
            https://en.wikipedia.org/wiki/Window_function#Planck-taper_window
        """
        # construct a Planck taper window
        w = np.ones(N)
        if nleft:
            w[0] *= 0
            zleft = np.array([nleft * (1./k + 1./(k-nleft))
                                for k in range(1, nleft)])
            w[1:nleft] *= expit(-zleft)
        if nright:
            w[N-1] *= 0
            zright = np.array([-nright * (1./(k-nright) + 1./k)
                                for k in range(1, nright)])
            w[N-nright:N-1] *= expit(-zright)
        return w

    return gwpy_planck(N, nleft=nleft, nright=nright)

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


