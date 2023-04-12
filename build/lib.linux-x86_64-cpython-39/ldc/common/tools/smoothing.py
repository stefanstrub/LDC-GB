import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d, splrep, splev

#pylint:disable=C0103

def moving_average(x, w):
    """ Smooth x signal with a window of soze w.
    """
    return np.convolve(x, np.ones(w), 'valid') / w


def mean_polyfit(PX, freq, Nmean, fmid=1.e-3):
    """ Smooth frequency domain PX signal using a loglog polynomial fit over the mean signal.
    """
    Smean = moving_average(PX, Nmean)
    fr_mean = freq[int(Nmean/2):-int(Nmean/2)]
    iem = len(fr_mean)
    if fr_mean[-1]>0.026:
        iem = np.argwhere(fr_mean>0.026)[0][0]
    ibm = np.argwhere(fr_mean>fmid)[0][0]
    z1 = np.polyfit(np.log10(fr_mean[:ibm]), np.log10(Smean[:ibm]), 30)
    z2 = np.polyfit(np.log10(fr_mean[ibm:iem]), np.log10(Smean[ibm:iem]), 30)
    Ssmooth1 = np.poly1d(z1)
    Ssmooth2 = np.poly1d(z2)
    return fr_mean[:iem], [Ssmooth1, Ssmooth2], Smean[:iem]

def sm_poly(spoly, freq, fmid=1.e-3):
    """
    """
    Sm = np.piecewise(freq,[freq<=fmid,freq>fmid],
                      [lambda x: 10**spoly[0](np.log10(x)),
                       lambda x: 10**spoly[1](np.log10(x))])
    return Sm

def sm_spline(Sspl, freq):
    """
    """
    Sm = 10**interpolate.splev(np.log10(freq), Sspl)
    return Sm

def mean_splinefit(PX, freq, Nmean, order):
    """ Smooth frequency domain PX signal using a loglog spline interpolator over the mean signal.
    """
    Smean = moving_average(PX, Nmean)
    fr_mean = freq[int(Nmean/2):-int(Nmean/2)]
    iem = len(fr_mean)
    if fr_mean[-1]>0.026:
        iem = np.argwhere(fr_mean>0.026)[0][0]
    x = fr_mean[:iem]
    Y = Smean[:iem]
    logx = np.log10(np.logspace(np.log10(x[0]), np.log10(x[-1]), 800))
    logy = np.interp(logx, np.log10(x), np.log10(Y))
    func = interpolate.interp1d(logx,logy)
    yi   = func(np.log10(x))
    tcl  = interpolate.splrep(np.log10(x), yi, s=order)
    return fr_mean[:iem], tcl, Smean[:iem]

def logpolyfit(x, y, order):
    """ Polynomial fit in loglog.
    """
    p = np.polyfit(np.log10(x), np.log10(y), order)
    Z = np.poly1d(p)
    return np.power(10, Z(np.log10(x)))

def logsplinefit(x, y, order):
    """ Spline interpolation in loglog.
    """
    logx = np.log10(np.logspace(np.log10(x[0]), np.log10(x[-1]), 200))
    logy = np.interp(logx, np.log10(x), np.log10(y))
    func = interp1d(logx,logy)
    yi   = func(np.log10(x))
    tcl  = splrep(np.log10(x), yi, s=order)
    ys   = splev(np.log10(x), tcl)
    return np.power(10, ys)
