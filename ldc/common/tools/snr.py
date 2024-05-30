""" Tools to compute SNR from tdi data and noise PSD
"""

#pylint:disable=C0103

import numpy as np
from ldc.lisa.noise import Noise

def compute_tdi_snr(source, noise, AET=False, data=None, fmin=-1, fmax=-1, **kwargs):
    """Compute SNR from TDI

    noise can be a Noise object returned by
    ldc.lisa.noise.get_noise_model or a dataset of PSD FrequencySeries

    source and data are datasets of TDI FrequencySeries
    """
    if data is None:
        data = source

    k1 = "X" if not AET else "A"

    if fmin == -1:
        fmin = source[k1].f[0] if source[k1].f[0]>0 else source[k1].f[1]
    if fmax == -1:
        fmax = source[k1].f[-1]

    if AET:
        return compute_tdi_snr_aet(source, noise, data,
                                   fmin, fmax, **kwargs)
    return compute_tdi_snr_xyz(source, noise, data,
                               fmin, fmax,**kwargs)

def compute_tdi_snr_aet(source, noise, data, fmin, fmax, ignore_T=False, tdi2=False,
                        full_output=False):
    """ Compute SNR from TDI A,E,T.
    """
    # inverse covariance matrix
    freq = np.array(source["A"].sel(f=slice(fmin, fmax)).f)
    df = source["A"].attrs["df"]
    if isinstance(noise, Noise):
        SA = noise.psd(freq=freq, option='A', tdi2=tdi2)
        if not ignore_T:
            ST = noise.psd(freq=freq, option='T', tdi2=tdi2)
    else:
        SA = noise["A"].sel(f=freq, method="nearest").values
        if not ignore_T:
            ST = noise["T"].sel(f=freq, method="nearest").values

    dsnr = dict()
    snr_tot = 0  # total SNR
    cumsum = np.zeros((len(freq)))

    keys, arrs = (["A", "E", "T"], [SA, SA, ST]) if not ignore_T else (["A", "E"], [SA, SA])
    
    for k, SN in zip(keys, arrs):
        d = data[k].sel(f=freq, method="nearest").values
        s = source[k].sel(f=freq, method="nearest").values

        snr = np.nansum(np.real(d*np.conj(s)/SN))
        snr *= 4.0*df
        dsnr[k+"2"] = snr

        snr_tot += snr
        if full_output:
            cumsum+= np.nancumsum(np.real(d * np.conj(s))/SN)

    dsnr["tot2"] = snr_tot
    if full_output:
        dsnr["cumsum"] = cumsum* 4.*df
        dsnr["freq"] = freq
    return dsnr


def compute_tdi_snr_xyz(source, noise, data, fmin, fmax, tdi2=False, full_output=False):
    """ Compute SNR from TDI X,Y,Z.
    """

    # inverse covariance matrix
    freq = np.array(source["X"].sel(f=slice(fmin, fmax)).f)
    df = source["X"].attrs["df"]
    if isinstance(noise, Noise):
        SXX = noise.psd(freq=freq, option='X', tdi2=tdi2)
        SXY = noise.psd(freq=freq, option='XY', tdi2=tdi2)
    else:
        SXX = noise["X"].sel(f=freq, method="nearest").values
        SXY = noise["XY"].sel(f=freq, method="nearest").values
    Efact = (SXX*SXX+SXX*SXY-2*SXY*SXY)
    Efact[Efact==0] = np.inf
    Efact = 1/Efact
    EXX = (SXX+SXY)*Efact
    EXY = -SXY*Efact

    dsnr = dict()
    snr_tot = 0  # total SNR
    cumsum = np.zeros((len(freq)))

    for k in ["X", "Y", "Z"]:     # individual SNR
        d = data[k].sel(f=freq, method="nearest").values
        s = source[k].sel(f=freq, method="nearest").values

        snr = np.nansum(np.real(d*np.conj(s)/SXX))
        snr *= 4.0*df
        dsnr[k+"2"] = snr

        snr = np.sum(np.real(d*np.conj(s)*EXX))
        snr_tot += snr
        if full_output:
            cumsum+= np.nancumsum(np.real(d * np.conj(s)*EXX))
    for k1,k2 in [("X", "Y"), ("X", "Z"), ("Y", "Z")]:
        d1 = data[k1].sel(f=freq, method="nearest").values
        s1 = source[k1].sel(f=freq, method="nearest").values
        d2 = data[k2].sel(f=freq, method="nearest").values
        s2 = source[k2].sel(f=freq, method="nearest").values

        snr = np.sum(np.real(d1*np.conj(s2)*EXY))
        snr+= np.sum(np.real(d2*np.conj(s1)*EXY))
        snr_tot += snr
        if full_output:
            cumsum+= np.cumsum(np.real(d1 * np.conj(s2)*EXY))
            cumsum+= np.cumsum(np.real(d2 * np.conj(s1)*EXY))
    snr_tot *= 4.*df
    dsnr["tot2"] = snr_tot

    if full_output:
        dsnr["cumsum"] = cumsum* 4.*df
        dsnr["freq"] = freq
        dsnr["EXX"] = EXX
        dsnr["EXY"] = EXY
    return dsnr
