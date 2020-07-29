import numpy as np

def get_freq_indices(f1, f2, fmin=-1, fmax=-1, i=0):
    """ Return indices of a common frequency range.  
    """

    ihmax = np.argwhere(f1 >=fmax)[0][0] if fmax>0 else len(f1)-1
    ihmin = np.argwhere(f1 >=fmin)[0][0] if fmin>0 else 1
    idmax = np.argwhere(f2 >=fmax)[0][0] if fmax>0 else len(f2)-1
    idmin = np.argwhere(f2 >=fmin)[0][0] if fmin>0 else 1

    if (f1[ihmin] != f2[idmin] or f1[ihmax] != f2[idmax]) and i==0:
        fmin = max(f1[ihmin], f2[idmin])
        fmax = min(f1[ihmax], f2[idmax])
        return get_freq_indices(f1, f2, fmin=fmin, fmax=fmax, i=1)
    else:
        return ihmin,ihmax,idmin,idmax

def compute_tdi_snr(source, noise, data=None, fmin=-1, fmax=-1, full_output=False):
    """ Compute SNR from TDI X,Y,Z. 
    
    noise is a Noise object, return by ldc.lisa.noise.get_noise_model
    source ans data are FrequencySeries
    """
    if data is None:
        data = source

    # freq indices selection
    freq = np.array(source["X"].f)
    df = source["X"].attrs["df"]
    ihmin,ihmax,idmin,idmax = get_freq_indices(freq, np.array(data["X"].f), fmin=fmin, fmax=fmax)
    
    # inverse covariance matrix
    SXX = noise.psd(freq=freq, option='X')
    SXY = noise.psd(freq=freq, option='XY')
    SXX, SXY = SXX[ihmin:ihmax], SXY[ihmin:ihmax]
    Efact = 1/(SXX*SXX+SXX*SXY-2*SXY*SXY)
    EXX = (SXX+SXY)*Efact
    EXY = -SXY*Efact

    dsnr = dict()
    for k in ["X", "Y", "Z"]:     # individual SNR
        snr = np.sum(np.real(data[k][idmin:idmax]*np.conj(source[k][ihmin:ihmax])/SXX))
        snr *= 4.0*df
        dsnr[k+"2"] = snr

    snr_tot = 0  # total SNR
    cumsum = np.zeros((len(freq)))
    for k in ["X", "Y", "Z"]:
        snr = np.sum(np.real(data[k][idmin:idmax]*np.conj(source[k][ihmin:ihmax])*EXX))
        snr_tot += snr
        if full_output:
            cumsum[ihmin:ihmax]+= np.cumsum(np.real(data[k][idmin:idmax] *\
                                                    np.conj(source[k][ihmin:ihmax])*EXX))
    for k1,k2 in [("X", "Y"), ("X", "Z"), ("Y", "Z")]:
        snr = np.sum(np.real(data[k1][idmin:idmax]*np.conj(source[k2][ihmin:ihmax])*EXY))
        snr+= np.sum(np.real(data[k2][idmin:idmax]*np.conj(source[k1][ihmin:ihmax])*EXY))
        snr_tot += snr
        if full_output:
            cumsum[ihmin:ihmax]+= np.cumsum(np.real(data[k1][idmin:idmax] *\
                                                    np.conj(source[k2][ihmin:ihmax])*EXY))
            cumsum[ihmin:ihmax]+= np.cumsum(np.real(data[k2][idmin:idmax] *\
                                                    np.conj(source[k1][ihmin:ihmax])*EXY))
    snr_tot *= 4.*df
    dsnr["tot2"] = snr_tot
    if full_output:
        dsnr["cumsum"] = cumsum* 4.*df
        dsnr["freq"] = freq[ihmin:ihmax]
        dsnr["EXX"] = EXX
        dsnr["EXY"] = EXY
    return dsnr

