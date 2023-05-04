"""Evaluation of galactic binaries submission for the sangria data
challenge.

Authors: M. Le Jeune lejeune@apc.in2p3.fr
"""
import copy
import glob
import numpy as np
import numpy.lib.recfunctions as recf
import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from tqdm import tqdm
import yaml

from data_tools import get_light_gb_injection, mbhb_free_data, download
from ldc.lisa.noise import AnalyticNoise
from ldc.common.series import XYZ2AET, TDI, FrequencySeries
from ldc.waveform.fastGB import FastGB
from ldc.common.tools import compute_tdi_snr
from ldc.common.tools import window

from lisacattools.catalog import GWCatalogs, GWCatalogType


dname = {'Frequency':'Frequency', 
         'Frequency Derivative':'FrequencyDerivative',
         'Amplitude':'Amplitude', 'Ecliptic Longitude':'EclipticLongitude',
         'Initial Phase':'InitialPhase', 'Polarization':'Polarization', 
         'Ecliptic Latitude':'EclipticLatitude',
         'Inclination':'Inclination'}

class GBEval:

    def __init__(self, team, workdir, sangria_fn='LDC2_sangria_training_v2.h5',
                 submitted_noise=False):
        """Load catalogs and noise model. 
        """
        self.workdir = workdir
        self.team = team
        self.teamdir = os.path.join(workdir, team)
        self.inj_cat = self.load_injection(sangria_fn=os.path.join(workdir, sangria_fn)) 
        self.inj_cat = self.inj_cat[self.inj_cat["snr"]>4]
        self.sub_cat = self.load_submission('catalog')
        self.noise, self.noise_name = self.load_noise(submitted_noise)
        tdidata = TDI.load(os.path.join(workdir, sangria_fn), "obs/tdi")
        dt = tdidata.X.attrs["dt"]
        df = 1/len(tdidata.t)/dt
        self.GB = FastGB(delta_t=dt, T=1/df)

    def get_pdf(self, i_inj=None, i_sub=None, full_chain=False):
        """Return source pdf as recarray, either by injection index or
        submisison index.

        """
        if i_inj is not None:
            if not self.match[i_inj]:
                print("no match")
                return
            i_c = np.argmax(np.array(self.match[i_inj])[:,1])
            i_sub = int(np.array(self.match[i_inj])[i_c,0])
        if self.team in ['apc-l2it', 'eth']:
            f0 = self.sub_cat["Frequency"][i_sub]
            idx = np.argmin(np.abs(self.pdf_f0-f0*1000))
            csv_file = self.pdf_file[idx]
            ldc_chain = np.genfromtxt(csv_file, names=True, delimiter=",")
            ldc_chain["Frequency"] *= 1e-3
            negative_longitude_mask = ldc_chain["EclipticLongitude"] < 0
            ldc_chain["EclipticLongitude"][negative_longitude_mask] += 2*np.pi
            return ldc_chain

        if self.team=='msfc-montana':
            record = self.sub_cat["record"][i_sub]
            ucb_samples_attr = self.ucb_catalog.get_attr_source_samples(record)
            ucb_samples = self.ucb_catalog.get_source_samples(record, ucb_samples_attr)
            pdf = ucb_samples.to_records(index=False)
            new_names = list(pdf.dtype.names).copy()
            for n in dname.keys():
                if n in new_names:
                    i = new_names.index(n)
                    new_names[i] = dname[n]
            pdf.dtype.names = new_names
            return pdf
            
    def load_submission(self, what):
        """ Convert yml into a recarray
        """
        if self.team=="eth":
            yml_solution = os.path.join(self.teamdir, "ETH-LDC2-sangria-training-v2-training-gb-only-SNR7.yaml")
            if not os.path.exists(yml_solution):
                download(yml_solution)
            csv_dir = os.path.join(self.teamdir, "Chain")
            if not os.path.exists(csv_dir):
                download(csv_dir)
            self.pdf_file = sorted(glob.glob(os.path.join(csv_dir, "*.csv")))
            self.pdf_f0 = [float(fn.split('frequency')[-1].split('nHz')[0])/10**6 for fn in self.pdf_file]
            if what=='catalog':
                L = yaml.load(open(yml_solution, "r"), yaml.Loader)
                df = pd.DataFrame(L["estimates"])
                names = list(df.columns)
                cat = np.rec.fromarrays([np.nan*np.zeros((len(df)))]*len(names),
                                        names=names)
                for i_n, n in enumerate(df.columns):
                    v = [df[n][i] for i in range(len(df))]
                    v = np.array([float(e) for e in v])
                    cat[n] = v[:]
                units = L["units"]
                if units["Frequency"]=='mHz':
                    cat["Frequency"] *= 1e-3
                negative_longitude_mask = cat["EclipticLongitude"] < 0
                cat["EclipticLongitude"][negative_longitude_mask] += 2*np.pi
                idx2 = np.arange(len(df))
                cat = recf.append_fields(cat, ['index'], [idx2], usemask=False)
                return cat

        if self.team=="apc-l2it":
            yml_solution = os.path.join(self.teamdir, "lejeune-ldc2a-training-gb-new.yaml")
            if not os.path.exists(yml_solution):
                download(yml_solution)
            csv_dir = os.path.join(self.teamdir, "pdf-new")
            if not os.path.exists(csv_dir):
                download(csv_dir)
            self.pdf_file = sorted(glob.glob(os.path.join(csv_dir, "*.csv")))
            pdf_fmin = [float(fn.split('-')[-2]) for fn in self.pdf_file]
            pdf_fmax = [float(fn.split('-')[-1].split(".pdf")[0]) for fn in self.pdf_file]
            self.pdf_f0 = (np.array(pdf_fmin)+np.array(pdf_fmax))/2.
            if what=='catalog':
                L = yaml.load(open(yml_solution, "r"), yaml.Loader)
                df = pd.DataFrame(L["estimates"])
                names_err = [n+'_error' for n in df.columns]
                names = list(df.columns)+names_err
                cat = np.rec.fromarrays([np.nan*np.zeros((len(df)))]*len(names),
                                        names=names)
                for i_n, n in enumerate(df.columns):
                    v = [df[n][i].split(" +/- ") for i in range(len(df))]
                    v = np.array([[float(e[0]), float(e[1])] for e in v])
                    cat[n] = v[:,0]
                    cat[names_err[i_n]] = v[:,1]
                units = L["units"]
                if units["Frequency"]=='mHz':
                    cat["Frequency"] *= 1e-3
                    cat["Frequency_error"] *= 1e-3
                idx2 = np.arange(len(df))
                cat = recf.append_fields(cat, ['index'], [idx2], usemask=False)
                return cat

        if self.team=="msfc-montana" and what=='catalog':
            ucb_fn = os.path.join(self.teamdir, 'UCB.h5')
            if not os.path.exists(ucb_fn):
                download(ucb_fn)
            catalogs = GWCatalogs.create(GWCatalogType.UCB, self.teamdir, "UCB.h5")
            self.ucb_catalog = catalogs.get_last_catalog()
            T = self.ucb_catalog.get_dataset("detections")
            #T = T[T["evidence"]>0.9]
            dname["evidence"] = 'evidence'
            N = len(dname.keys())
            cat = np.rec.fromarrays([np.nan*np.zeros((len(T)))]*N,
                                    names=list(dname.values()))
            for k,v in dname.items():
                cat[v] = T[k]
            idx1 = [i for i in T.index]
            idx2 = np.arange(len(T))
            cat = recf.append_fields(cat, ['record', 'index'],
                                     [np.array(idx1), idx2], usemask=False)
            return cat
        

    def load_noise(self, submitted_noise=False):
        """ Return noise model as a psd_a(f). 
        """
        if submitted_noise:
            if self.team=='eth_radler_noise':
                noise_model = "SciRDv1"
                frequencies = np.logspace(-4., -1.9, 1000)
                ldc_noise = AnalyticNoise(frequencies, model=noise_model, wd=1)
                f = spline(frequencies, ldc_noise.psd(option='A'))
                return f, noise_model
            if self.team=='eth':
                noise_fn = os.path.join(self.teamdir, 'ETH_sangria_noise.csv')
                if not os.path.exists(noise_fn):
                    download(noise_fn)
                psd = pd.read_csv(noise_fn, delimiter=",")  
                f = spline(psd['f'], psd['A'])
                return f, self.team
            if self.team=='apc-l2it':
                noise_fn = os.path.join(self.teamdir, 'APC_L2IT_sangria-noise.csv')
                if not os.path.exists(noise_fn):
                    download(noise_fn)
                psd = np.genfromtxt(noise_fn, names=True, delimiter=",")['meanPSD']
                frqs = np.logspace(-4., -1.9, 100)
                f = spline(frqs, psd)
                return f, self.team
            if self.team=='msfc-montana':
                noise_fn = os.path.join(self.teamdir, 'psd.dat')
                psd = np.genfromtxt(noise_fn, names=True)
                f = spline(psd['f'][0:-1:4], 2*psd["SnA2"][0:-1:4])
                return f, self.team
                
        frequencies = np.logspace(-4., -1.9, 1000)
        ldc_noise = AnalyticNoise(frequencies, model="sangria", wd=1)
        f = spline(frequencies, ldc_noise.psd(option='A'))
        return f, 'sangria'

    def load_injection(self, sangria_fn, pop=['dgb', 'igb', 'vgb'], snr_thresh=2):
        """Load a light version of the true catalog used to build the dataset,
        based on some SNR threshold. .
        """
        if not os.path.exists(sangria_fn):
            download(sangria_fn)
        pop_str = '-'.join(pop)
        light_injection = os.path.join(self.workdir,
                                       f"gb_light_catalog.npy")
        if not os.path.exists(light_injection):
            cat = get_light_gb_injection(sangria_fn)
            np.save(light_injection, cat)
        else:
            cat = np.load(light_injection)
        for p in ["dgb", "igb", "vgb"]:
            if p not in pop:
                cat = cat[~cat["type"]==p]
        return cat[cat["snr"]>snr_thresh]
        

    def compute_overlap(self, s1, s2, oversample=2, A1=None, E1=None, returnAE=False):
        """Compute overlap between 2 sources, weighted by noise.

        """
        if A1 is None:
            A1,E1,T1 = XYZ2AET(*self.GB.get_fd_tdixyz(template=s1, oversample=oversample))
        A2,E2,T2 = XYZ2AET(*self.GB.get_fd_tdixyz(template=s2, oversample=oversample))
        SA = np.mean(self.noise(A1.f))

        # use explicit indices to avoid precision issues with f
        size1, size2 = len(A1), len(A2)
        kmin1, kmin2 = A1.attrs["kmin"], A2.attrs['kmin']

        bd, bm = (0, kmin1-kmin2) if kmin2<kmin1 else (kmin2-kmin1, 0)
        ed = size1 if size1-bd<size2-bm else size2-bm+bd
        em = bm+(ed-bd)

        olap = 0
        if ed>bd and em>bm:
            n1 = np.sum((np.abs(A1.values)**2 + np.abs(E1.values)**2) /SA)
            n2 = np.sum((np.abs(A2.values)**2 + np.abs(E2.values)**2) /SA)
            AAi = A1.values[bd:ed]*np.conj(A2.values[bm:em]) # reduced to common freq
            EEi = E1.values[bd:ed]*np.conj(E2.values[bm:em])
            olap = np.sum(np.real(AAi) + np.real(EEi))/SA
            olap = np.abs(olap)/np.sqrt(n1*n2)
        if returnAE:
            return olap, A1, E1, A2, E2
        return olap
        
    def make_match(self, thresh=0.1, oversample=2):
        """Match submitted catalog with the reference one. 
        
        For each source in the reference catalog, find one or several
        match(es) in the submitted catalog (limit to f+/- ~10uHz), for
        a given correlation threshold. 
        
        keys:injection index
        values:[(submission index, overlap)]
        """
        add_field = ["snr"]
        match = {}
        deltaf = 1e-5
        for i_s in tqdm(range(len(self.inj_cat))):
            s = self.inj_cat[i_s]
            pGB1 = dict(zip(s.dtype.names, s))
            match[i_s] = []
            f0 = s["Frequency"]
            selec = (self.sub_cat["Frequency"]>f0-deltaf)
            selec &= (self.sub_cat["Frequency"]<f0+deltaf)
            A1,E1,T1 = XYZ2AET(*self.GB.get_fd_tdixyz(template=pGB1,
                                                      oversample=oversample))
            for s1 in self.sub_cat[selec]:
                pGB2 = dict(zip(s1.dtype.names, s1))
                olap = self.compute_overlap(pGB1, pGB2, A1=A1, E1=E1,
                                            oversample=oversample)
                if olap>thresh:
                    match[i_s].append((s1["index"], olap))
        return match

    def get_AE(self, s, oversample=2):
        """ Return TDI A and E for a given source. 
        """
        pGB = dict(zip(s.dtype.names, s))
        X,Y,Z = self.GB.get_fd_tdixyz(template=pGB, oversample=oversample)
        A,E,T = XYZ2AET(X, Y, Z)
        A.attrs = X.attrs
        return A, E
    
    def compute_snr(self, s, oversample=2):
        """ <s|s>
        """
        pGB = dict(zip(s.dtype.names, s))
        X,Y,Z = self.GB.get_fd_tdixyz(template=pGB, oversample=oversample)
        A,E,T = XYZ2AET(X, Y, Z)
        A.attrs = X.attrs
        noise = {"A":FrequencySeries(self.noise(X.f), fs=X.f)}
        snr = compute_tdi_snr(dict({"A":A,"E":E,"T":T}), noise, AET=True,
                              ignore_T=True)
        snr = np.sqrt(snr['A2']+snr['E2'])
        return snr

    
    def add_submission_snr(self):
        """Add snr column to submitted catalog, computed from given noise
        model.
        """
        snrl = np.zeros((len(self.sub_cat)))
        for  i in tqdm(range(len(self.sub_cat))):
            snrl[i] = self.compute_snr(self.sub_cat[i])
        self.sub_cat = recf.append_fields(self.sub_cat, ["snr"], [snrl], usemask=False)

    def get_multimatch(self, thresh=0.9):
        """Return multi-match, ie where an injection match several
        sources in the submitted results.
        
        Only overlap above thresh are considered as a match here. 
        Submitted sources which have another match are discarded. 
        """
        multim_inj = [k for k,v in self.match.items() if len(v)>1]
        multim_inj = [k for k in multim_inj
                      if (np.array(self.match[k])[:,1]>thresh).sum()>1]

        discard = np.zeros((len(multim_inj))).astype(bool)
        for k, inj in enumerate(multim_inj):
            f0 = self.inj_cat["Frequency"][inj]
            neighbours = np.where((self.inj_cat["Frequency"]>f0-1e-5) & (self.inj_cat["Frequency"]<f0+1e-5))[0]
            neighbours = list(neighbours)
            neighbours.remove(inj)
            neighbours_match = [j[0] for i in neighbours
                                for j in self.match[i] if j[1]>thresh]

            detections = np.array(self.match[inj])[:,1]>thresh
            detections = np.array(self.match[inj])[detections, 0]
            for n in detections:
                if n in neighbours_match:
                    discard[k] = True
        multim_inj = np.array(multim_inj)[~discard]
        return multim_inj

    
    def get_lowmatch(self, thresh=0.9):
        """Return low-match, ie where a recovered source fits a
        combination of injections. 
        """
        sub, inj = self._get_best_match(thresh=thresh)
        low_match = {}
        for s in sub:
            wh = np.where(sub["index"]==s["index"])[0]
            if len(wh)>1:
                low_match[int(s["index"])] = inj[wh]
        return low_match

    def get_spurious_detection(self):
        """Return submitted sources which don't match any
        injection (with correlation below 0.1).
        """
        all_match = list(set([m[0] for k,v in self.match.items()
                              for m in v]))
        false = [i for i in self.sub_cat["index"] if i not in all_match]
        false = np.array(false).astype(int)
        return self.sub_cat[false]

    def get_partial_detection(self, thresh=0.9, injection=True):
        """Return submitted sources which partially match an
        injection (with correlation below thresh).
        """
        if injection:
            partial = [k for k,v in self.match.items()
                       if len(v)>0 and np.array(v)[:,1].max()<thresh]
            return self.inj_cat[partial]
        all_match = list(set([m[0] for k,v in self.match.items()
                              for m in v if m[1]>thresh]))
        no_no_match = list(set([m[0] for k,v in self.match.items()
                                for m in v]))
        false = [i for i in self.sub_cat["index"]
                 if i not in all_match and i in no_no_match]
        false = np.array(false).astype(int)
        return self.sub_cat[false]

    def get_max_overlap(self):
        """Return max overlap found for each injection.
        """
        max_olap = [np.array(v)[:,1].max() if len(v)>0 else 0 for k,v in self.match.items()]
        return max_olap
        
    def get_true_detection(self, thresh=0.9, injection=True):
        """Return submitted sources which have at least one match with an
        injection (with correlation above thresh).
        """
        if injection:
            all_match = [k for k,v in self.match.items()
                         if len(v)>0 and np.array(v)[:,1].max()>thresh]
            return self.inj_cat[all_match]
        all_match_ = list(set([m[0] for k,v in self.match.items()
                               for m in v if m[1]>thresh]))
        return self.sub_cat[np.array(all_match_).astype(int)]

    def _get_best_match(self, thresh=0.9, return_index=False):
        """Find best match and return aligned catalogs.
        """
        all_match = [np.array(v)[np.argmax(np.array(v)[:,1]),0]
                     for k,v in self.match.items()
                     if len(v)>0 and np.array(v)[:,1].max()>thresh]
        all_match = np.array(all_match).astype(int)
        all_ref = [k for k,v in self.match.items()
                   if len(v)>0 and np.array(v)[:,1].max()>thresh]
        sub, inj = self.sub_cat[all_match], self.inj_cat[all_ref]
        if return_index:
            return sub, inj, all_ref
        return sub, inj
    
    def get_single_detection(self, thresh=0.9):
        """Find one to one match and return submitted and corresponding
        injected sources
        """
        sub, inj = self._get_best_match(thresh=thresh)
        inj, idx = np.unique(inj, return_index=True)
        sub = sub[idx]
        sub, idx = np.unique(sub, return_index=True)
        inj = inj[idx]
        return inj, sub
    
    def get_missing(self, thresh=0.1):
        """Return injections which don't have been recovered at all
        (with correlation below thresh).
        """
        missing = [k for k,v in self.match.items()
                   if len(v)==0 or np.array(v)[:,1].max()<thresh]
        return self.inj_cat[missing]
    
    def save_workspace(self):
        """ Save all matched quantities. 
        """
        d = dict({"noise":self.noise,
                  "inj_cat":self.inj_cat, "sub_cat":self.sub_cat})
        if hasattr(self, 'match'):
                 d["match"] = self.match
        fn = os.path.join(self.teamdir, f"match_noise-{self.noise_name}.pkl")
        pickle.dump(d, open(fn, "wb"))

    def load_from_workspace(self):
        """ Load existing match or recompute it. 
        """
        fn = os.path.join(self.teamdir, f"match_noise-{self.noise_name}.pkl")
        if os.path.exists(fn):
            d = pickle.load(open(fn, "rb"))
            self.match = d["match"]
            self.sub_cat = d["sub_cat"]
            self.inj_cat = d["inj_cat"]
        else:
            match = self.make_match()
            self.match = match
            self.add_submission_snr()
            self.save_workspace() 

if __name__ == '__main__':        

    workdir = "/home/stefan/LDC/Sangria/evaluation"
    # gb_apc = GBEval('apc-l2it', workdir, submitted_noise=True)
    # gb_usa = GBEval('msfc-montana', workdir, submitted_noise=True)
    gb_eth = GBEval('eth', workdir, submitted_noise=True)

    # gb_apc.load_from_workspace()
    # gb_usa.load_from_workspace()

    gb_eth.load_from_workspace()

    if 0:
        samples_list = list()
        for i_sub in range(len(gb_apc.sub_cat)):
            pdf = gb_apc.get_pdf(i_sub=i_sub)
            if pdf is not None:
                names = list(pdf.dtype.names).copy()
                names[names.index("EclipticLatitude")] = "Ecliptic Latitude"
                names[names.index("EclipticLongitude")] = "Ecliptic Longitude"
                pdf.dtype.names = names
                pdf = pd.DataFrame(pdf)
                samples_list.append(pdf[["Ecliptic Latitude", "Ecliptic Longitude"]])
            else:
                print('missing', i_sub)

        all_sources = pd.concat(samples_list)
        from lisacattools import HPhist
        nside = 64
        hpmap = HPhist(all_sources, nside)
        fig = plt.figure(figsize=(8, 6), dpi=100)
        ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection="geo degrees mollweide")
        ax.grid()
        # use logarithmic scaling for density
        ax.imshow_hpx(np.log10(hpmap + 1), cmap="plasma")

    max_olap = np.array(gb_usa.get_max_overlap())
    issue = (max_olap==0) & (gb_usa.inj_cat["snr"]>100)
    idx = np.where(issue)
    print(gb_usa.inj_cat[issue])

  
