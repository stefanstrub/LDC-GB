"""Evaluation of galactic binaries submission for the sangria data
challenge.

Authors: M. Le Jeune lejeune@apc.in2p3.fr
"""

import numpy as np
import numpy.lib.recfunctions as recf

from ldc.lisa.noise import AnalyticNoise, simple_snr
from ldc.common.series import XYZ2AET, TDI
from ldc.common.tools import compute_tdi_snr
from ldc.waveform.fastGB import FastGB
import ldc.io.hdf5 as h5io
from ldc.common.tools import window

from ldc.waveform.lisabeta import FastBHB
"""
TODO
- check simple snr
"""
import os
import yaml

def download(fn):
    """ Download from ldc website.
    """
    dirname = os.path.dirname(fn)
    basename = os.path.basename(fn)

    url = dict({'LDC2_sangria_training_v2.h5':'https://lisa-ldc.lal.in2p3.fr/media/uploads/LDC2_sangria_training_v2.h5',
                'lejeune-ldc2a-training-gb-new.yaml':'https://lisa-ldc.lal.in2p3.fr/media/uploads/submission/lejeune-ldc2a-training-gb-new.yaml',
                'pdf-new':'https://lisa-ldc.lal.in2p3.fr/media/uploads/submission/pdf-new.tar.gz',
                'UCB.h5':'https://lisa-ldc.lal.in2p3.fr/media/uploads/submission/MSFC-Montana_12mo_catalog.tgz',
                'APC_L2IT_sangria-noise.csv':'https://lisa-ldc.lal.in2p3.fr/media/uploads/submission/APC_L2IT_noise_MBHBs.tar.gz',
                'APC_L2IT-ldc2a-mbhb.yaml':'https://lisa-ldc.lal.in2p3.fr/media/uploads/submission/APC_L2IT_noise_MBHBs.tar.gz',
                })
    o = os.path.join(dirname, url[basename].split('/')[-1])
    os.system(f'wget {url[basename]} -O {o}')
    if o[-2:]=='gz':
        os.system(f"tar -xzf {o} -C {fn}")


def mbhb_data(sangria_fn, team, workdir):
    """ Return MBHB tdi A,E 
    """
    tdidata = TDI.load(sangria_fn, "obs/tdi")
    dt = tdidata.X.attrs["dt"]
    T = len(tdidata.t)
    freqs = np.fft.rfftfreq(T, d=dt)
    A_data = np.zeros(len(freqs), dtype=np.complex128)
    E_data = np.zeros(len(freqs), dtype=np.complex128)
    T_data = np.zeros(len(freqs), dtype=np.complex128)

    teamdir = os.path.join(workdir, team)
    if team=="apc-l2it":
        mbhb_param = os.path.join(teamdir, 'APC_L2IT-ldc2a-mbhb.yaml')
        L = yaml.load(open(mbhb_param, "r"), yaml.Loader)
        phenomD = FastBHB(approx="IMRPhenomD", T=T, delta_t=dt, bbh_type='mbhb')
        for k in L.keys():
            if k[0:3]=='src':
                param = L[k][0]
                param["Spin1"] = param["ProjectedSpin1"]
                param["Spin2"] = param["ProjectedSpin2"]
                A, E, T = phenomD.get_fd_tdiaet(template=param)
                kmin = A.attrs["kmin"]
                A_data[kmin:kmin+len(A.values)]+= A.values
                E_data[kmin:kmin+len(A.values)]+= E.values
                T_data[kmin:kmin+len(A.values)]+= T.values
    return A_data, E_data, T_data

def mbhb_free_data(sangria_fn):
    """Return sangria data in frequency domain, AET, with perfect removal
    of MBHB. 

    """
    tdidata = TDI.load(sangria_fn, "obs/tdi")
    attrs = tdidata.X.attrs
    mbhb = TDI.load(sangria_fn, "sky/mbhb/tdi")
    tdidata -= mbhb
    tdidata.X.attrs = attrs; tdidata.Y.attrs = attrs; tdidata.Z.attrs = attrs; 
    tdidata.XYZ2AET()
    tdidata = TDI(dict(zip(tdidata.keys(), 
                           [tdidata[k].ts.fft(win=window, margin=1e6, kap=4.e-6)
                            for k in tdidata.keys()])))
    return tdidata
    

def get_light_gb_injection(sangria_fn, snr_thresh=2):
    """Return a light version of the GB injection catalog, all source
    types in a single recarray.

    """
    pop = ['dgb', 'igb', 'vgb']
    p = ['Amplitude', 'EclipticLatitude', 'EclipticLongitude', 'Frequency', 
         'FrequencyDerivative',  'Inclination',  'InitialPhase', 'Polarization']
    cat = []
    stype = []
    for i,k in enumerate(pop):
        scat, units = h5io.load_array(sangria_fn, f"sky/{k}/cat")
        snr = simple_snr(scat['Frequency'], scat["Amplitude"], incl=scat["Inclination"])
        selec = snr>snr_thresh
        scat = scat[selec]
        
        filtered = [x for x in scat.dtype.names if x in p]
        scat = scat[filtered]
        if i>0:
            scat = scat.astype(cat[0].dtype)
        stype += [k]*len(scat)
        cat.append(scat)

    cat = np.hstack(cat)
    cat = recf.append_fields(cat, ["type"], [np.array(stype)], usemask=False)

    tdidata = TDI.load(sangria_fn, "obs/tdi")
    dt = tdidata.X.attrs["dt"]
    df = 1/len(tdidata.t)/dt
    frequencies = np.fft.rfftfreq( int(len(tdidata.t) / dt), d=dt)
    noise = AnalyticNoise(frequencies, model="sangria", wd=1)
    GB = FastGB(delta_t=dt, T=1/df)
        
    def compute_snr(s):
        pGB = dict(zip(s.dtype.names, s))
        X,Y,Z = GB.get_fd_tdixyz(template=pGB, oversample=2)
        A,E,T = XYZ2AET(X, Y, Z)
        A.attrs = X.attrs
        snr = compute_tdi_snr(dict({"A":A,"E":E,"T":T}), noise, AET=True, ignore_T=True)
        snr = np.sqrt(snr['A2']+snr['E2'])
        return snr
        
    snrl = [compute_snr(s) for s in cat]
    cat = recf.append_fields(cat, ["snr"], [snrl], usemask=False)
    return cat

