import numpy as np

def compute_tdi_snr(source, noise, data=None, fmin=-1, fmax=-1, full_output=False):
    """ Compute SNR from TDI X,Y,Z. 
    
    noise is a Noise object, return by ldc.lisa.noise.get_noise_model
    source ans data are FrequencySeries
    """
    if data is None:
        data = source

    # freq indices selection
    df = source["X"].attrs["df"]
    if fmin == -1:
        fmin = source["X"].f[0] if source["X"].f[0]>0 else source["X"].f[1] 
    if fmax == -1:
        fmax = source["X"].f[-1]
    
    # inverse covariance matrix
    freq = np.array(source["X"].sel(f=slice(fmin, fmax)).f)
    SXX = noise.psd(freq=freq, option='X')
    SXY = noise.psd(freq=freq, option='XY')
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

        snr = np.sum(np.real(d*np.conj(s)/SXX))
        snr *= 4.0*df
        dsnr[k+"2"] = snr

        snr = np.sum(np.real(d*np.conj(s)*EXX))
        snr_tot += snr
        if full_output:
            cumsum+= np.cumsum(np.real(d * np.conj(s)*EXX))
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

