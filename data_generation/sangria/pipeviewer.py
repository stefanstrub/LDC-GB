"""
"""
import h5py
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import os
from ldc.lisa.projection import from_file, ProjectedStrain
from ldc.lisa.noise import simple_snr
from ldc.common.series import TimeSeries, FrequencySeries, TDI
import ldc.io.hdf5 as h5io
from ldc.lisa.orbits import Orbits
from ldc.common.tools import window
import ldc.waveform.fastGB as fastGB
from ldc.waveform.waveform import HpHc
from ldc.waveform.lisabeta import FastBHB
from ldc.lisa.noise import get_noise_model

def get_psd(tdi, k='X'):
     win = window(np.array(tdi.t-tdi.t[0]))#, xl=5000)
     d = tdi[k].data
     d[np.isnan(d)] = 0
     dt = tdi.t.data[1]-tdi.t.data[0]
     f, psdX =  scipy.signal.welch(win*d, fs=1.0/dt, window='hanning', nperseg=256*256)
     return f, psdX

class PipeViewer:
    """Plot all the data produced by the data generation pipeline.
    """

    def __init__(self, path, subpath='run1'):
        self.path = path
        if not isinstance(subpath, dict):
            self.subpath = dict({'dgb':subpath, 'igb':subpath})
        else:
            self.subpath = subpath
        self.tdi2 = False
            
    @classmethod
    def spritz(cls, path, subpath=None):
        return SpritzViewer(path)
        
    def get_catalog(self, source):
        if source in ["igb", "dgb"]:
            cfile = os.path.join(self.path, f"../{source}-pipeline/{self.subpath[source]}",
                                 f'{source}.h5')
        else:
            cfile = os.path.join(self.path, f'{source}.h5')
            if not os.path.exists(cfile):
                cfile = os.path.join(self.path, f'{source}-cat.h5')
        cat, descr = h5io.load_array(cfile)
        return cat, descr
    
    def plot_catalog(self, source):
        """ Catalog histograms
        """
        cat, descr = self.get_catalog(source)

        if source in ['igb', 'dgb']:
            plt.figure()
            plt.hist(np.log10(cat["Amplitude"]), bins=100)
            plt.xlabel("log10 ampl")
            plt.yscale("log")
            plt.figure()
            plt.hist(cat["Frequency"], bins=100)
            plt.yscale("log"); plt.xscale("log");
            plt.xlabel("freq [Hz]")
            plt.figure()
            plt.hist(np.log(np.abs(cat["FrequencyDerivative"])), bins=100)
            plt.xlabel("log abs fdot")
            for k in ["Inclination", "EclipticLatitude", "EclipticLongitude", "InitialPhase", "Polarization"]:
                plt.figure(); plt.hist(cat[k], bins=100); plt.xlabel(k)
        if 'gb' in source:
            snr = simple_snr(cat['Frequency'], cat["Amplitude"], incl=cat["Inclination"])
            if source=="vgb":
                fig = plt.figure()
                ax = fig.add_subplot(111)
                plt.semilogy(snr, marker="o", ls='None')
                plt.ylabel("SNR")
                plt.xlabel("source index")
                for i_s, s in enumerate(cat):
                    ax.annotate(s["Name"], xy=(1.01*i_s, 1.01*snr[i_s]), textcoords='data', fontsize=8, rotation=45)
            else:
                plt.figure()
                plt.hist(snr, bins=100)
                plt.yscale("log");plt.xscale("log");
                plt.xlabel("snr")
        return cat

    def get_orbits(self, analytic=False):
        if not analytic:
            ofile = os.path.join(self.path, 'orbits.h5')
            orbits = Orbits.type(dict({'orbit_type':'file', 'filename':ofile, 'nominal_arm_length':2.5e9}))
        else:
            orbits = Orbits.type(dict({'orbit_type':'analytic', 'nominal_arm_length':2.5e9,
                                       "initial_position": 0, "initial_rotation": 0}))
        return orbits

    def get_config(self):
        return h5io.load_config(os.path.join(self.path, f"training.h5"), "sky/config")
   
    def plot_orbits(self):
        """ Plot orbits
        """
        orbits = get_orbits()
        p1 = orbits.compute_position(1, orbits.t)
        plt.figure()
        plt.plot(orbits.t, p1.T)
        plt.xlabel("Time [s]")
        plt.ylabel("Position of sc 1")
        print(f'tmax is {orbits.t[-1]}')
        tt = orbits.compute_travel_time(1, 2, orbits.t)
        plt.figure()
        plt.plot(orbits.t, tt)
        plt.xlabel("Time [s]")
        plt.ylabel("Travel time 1 to 2")
        return p1
                                     
        
    def load_strain(self, fn, link='1-2'):
        """ Load projected strains
        """
        yArm, source_names, links, t_min, t_max, dt = from_file(fn)
        ilink = links.index(link)
        yArm = TimeSeries(yArm[:,ilink], dt=dt, t0=t_min)
        return yArm
        
    def plot_strain(self, source, link='1-2', raw=False, time_domain=True, **kwargs):
        """ Plot projected strains
        """
        if source in ["igb", "dgb"] and raw:
            yfile = os.path.join(self.path, f"../{source}-pipeline/{self.subpath[source]}",
                                 f'{source}-y.h5')
        else:
            yfile = os.path.join(self.path, f'{source}-y.h5')

        if 'label' not in kwargs:
            label = source+'-'+link if not raw else 'raw '+source+'-'+link
            kwargs["label"] = label
        yArm = self.load_strain(yfile, link=link)
        if not time_domain:
            yArm = yArm.ts.fft()
            plt.loglog(yArm.f, np.abs(yArm), **kwargs)
            plt.xlabel("Freq [Hz]")
            plt.ylabel(f"|y| {link}")

        else:
            yArm.plot(**kwargs)
            plt.ylabel(f"y {link}")

    def get_tdi_x(self, source, raw=False, time_domain=True):
        """
        """
        if raw:
            tdifile = os.path.join(self.path, f'{source}-tdi-raw.h5')
            with h5py.File(tdifile, "r") as f:
                X = f['X'][:]
            ineg = X[:,0]>=0
            t = X[ineg,0]
            tdi = TimeSeries(X[ineg,1], t0=t[0], dt=t[1]-t[0])
        else:
            tdifile = os.path.join(self.path, f'{source}-tdi.h5')
            tdi, descr = h5io.load_array(tdifile)
            dt = tdi["t"][1]-tdi["t"][0]
            tdi = TimeSeries(tdi["X"], t0=tdi["t"][0], dt=dt)
            
        if not time_domain:
            tdi = tdi.ts.fft(win=window)

        return tdi

    def get_tt_order(self):
        #config = self.get_config()
        return 1#config["travel_time_order"]
   
    def get_tdi_x_from_y(self, source, t):
        #yfile = os.path.join(self.path, f"../{source}-pipeline/{self.subpath}",
        #                     f'{source}-y.h5')
        yfile = os.path.join(self.path, f'{source}-y.h5')
        p = ProjectedStrain(self.get_orbits())
        p.from_file(yfile)
        tdi_x = p.compute_tdi_x(t, tdi2=self.tdi2, tt_order=self.tt_order())
        return TimeSeries(tdi_x, t0=t[0], dt=t[1]-t[0])

    def get_tdi_x_old(self, source, time_domain=True):
        tdifile = os.path.join(self.path, f'{source}-lisanode/{source}-tdi.h5')
        tdi, descr = h5io.load_array(tdifile)
        t = tdi['X'][:,0]
        selec = t>=0
        t = t[selec]
        dt = t[1]-t[0]
        tdi_sangria_x = TimeSeries(tdi["X"][selec,1], dt=dt, t0=t[0])
        if not time_domain:
            tdi_sangria_x = tdi_sangria_x.ts.fft(win=window)
        return tdi_sangria_x
    
    def plot_tdi_x(self, source, time_domain=True, raw=False):
        """
        """
        tdi = self.get_tdi_x(source, raw=raw)
        if raw:
            tdi = tdi[0:60*60*24*30*4] # 1 month
        expected_tdi = self.get_tdi_x_from_y(source, tdi.t.values)
        win = window(tdi.t.values)
        tdi.data *= win
        expected_tdi *= win
        
        if time_domain:
            plt.subplot(211)
            tdi.plot(label=source)
            expected_tdi.plot(label="from y", color='k', ls='--', alpha=0.5)
            plt.ylabel("windowed TDI X")
            plt.legend(loc="upper right")
            plt.subplot(212)
            plt.plot(tdi.t, tdi.data-expected_tdi.data)
            plt.ylabel("Difference")
            plt.xlabel("Time [s]")
        else:
            tdi = tdi.ts.fft()
            expected_tdi = expected_tdi.ts.fft()
            plt.subplot(211)
            plt.semilogy(tdi.f, np.abs(tdi), label=source)
            plt.plot(expected_tdi.f, np.abs(expected_tdi), label="from y", color='k', ls='--', alpha=0.5)
            plt.ylabel("windowed TDI X")
            plt.subplot(212)
            plt.semilogy(tdi.f, np.abs(tdi.data-expected_tdi.data))
            plt.xlabel("Freq [Hz]")
            plt.ylabel("Difference")
        return tdi

    def check_mbhb_source(self, pMBHB=None, index=0, tdi=None, single=False):
        """ Check a particular MBHB source with lisabeta
        """
        if pMBHB is None:
            cfile = os.path.join(self.path, f'mbhb.h5')
            cat, descr = h5io.load_array(cfile)
            pMBHB = dict(zip(cat.dtype.names, cat[index]))
        if tdi is None:
            tdi = self.get_tdi_x('mbhb')
        dt = tdi.attrs['dt']
        tvec = tdi.t.data
        t_max = tvec[-1]+dt
        t_min = tvec[0]

        FBH = FastBHB("MBHB", T=t_max, delta_t=dt, approx='IMRPhenomD', orbits=self.get_orbits(analytic=True))
        #FBH.wvf_pars = waveform_params
        A,E,T = FBH.get_td_tdiaet(template=pMBHB)
        XYZ = TDI(dict(zip(["A", "E", "T"], [A, E, T])))
        XYZ.AET2XYZ()
        prange = slice(pMBHB["CoalescenceTime"]-1000, pMBHB["CoalescenceTime"]+500)
        sets = [XYZ.X.sel(t=prange), tdi.sel(t=prange)]
        labs = ['lisabeta', 'data']
        
        if single:
            p = ProjectedStrain(self.get_orbits(analytic=True))
            GW = HpHc.type("debug", "MBHB", "IMRPhenomD")
            GW.set_param(pMBHB)
            yArm = p.arm_response(t_min, t_max, dt, [GW], tt_order=self.get_tt_order())
            tdi_x = TimeSeries(p.compute_tdi_x(tvec, tdi2=self.tdi2, tt_order=self.get_tt_order()), dt=dt, t0=t_min)
            sets.append(tdi_x.sel(t=prange))
            labs.append('single')
            
        plt.subplot(2,1,1)
        for ts, l in zip(sets, labs):
            plt.plot(ts.t, ts, label=l)
        #plt.axis([pMBHB["CoalescenceTime"]-500, pMBHB["CoalescenceTime"]+500, None, None])
        plt.legend()

        plt.subplot(2,1,2)
        plt.semilogy(sets[0].t, np.abs(sets[0].data-sets[1].data), label=f'{labs[0]}-{labs[1]}')
        if len(sets)>2:
            plt.plot(sets[1].t, np.abs(sets[1].data-sets[2].data), label=f'{labs[1]}-{labs[2]}')
        plt.legend()
        #plt.axis([pMBHB["CoalescenceTime"]-500, pMBHB["CoalescenceTime"]+500, None, None])
        plt.grid()
        return XYZ

    def check_gb_source(self, pGB, source='dgb', tdi=None, single=False, pure_python=False):
        """ Check a particular GB source with fastGB
        """
        if tdi is None:
            tdi = self.get_tdi_x(source)
        dt = tdi.attrs['dt']
        t_min = tdi.attrs['t0']
        tvec = tdi.t.data
        t_max = tvec[-1]+dt
        tdi = tdi.ts.fft(win=window)
        GB = fastGB.FastGB(delta_t=dt, T=t_max, orbits=self.get_orbits(analytic=True)) # in seconds
        Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, radler=not(pure_python))
        tdi_fg = TDI(dict(zip(["X", "Y", "Z"], [Xs, Ys, Zs])))
        tdi_sub = tdi.sel(f=tdi_fg.f, method="nearest")

        sets = [[tdi_fg.X.real, tdi_sub.real],
                [tdi_fg.X.imag, tdi_sub.imag],
                [np.abs(tdi_fg.X), np.abs(tdi_sub)]]
        labs = ['fastGB', 'data']
        ylabs = ["X.real", "X.imag", "|X|"]

        if single:
            p = ProjectedStrain(self.get_orbits())
            GW = HpHc.type("debug", "GB", "TD_fdot")
            GW.set_param(pGB)
            yArm = p.arm_response(t_min, t_max, dt, [GW], tt_order=self.get_tt_order())
            tdi_x = TimeSeries(p.compute_tdi_x(tvec, tdi2=self.tdi2, tt_order=self.get_tt_order()), dt=dt, t0=t_min)
            tdi_x = tdi_x.ts.fft(win=window)
            tdi_x = tdi_x.sel(f=tdi_fg.f, method="nearest")
            
            sets[0].append(tdi_x.real)
            sets[1].append(tdi_x.imag)
            sets[2].append(np.abs(tdi_x))
            labs.append('single')
            

        i = 1
        for ts in sets:
            plt.subplot(2,3,i)
            for t1, l in zip(ts, labs):
                plt.plot(t1.f, t1, label=l, alpha=0.5)
            plt.axis([pGB["Frequency"]-1.5e-6, pGB["Frequency"]+1.5e-6, None, None])
            plt.legend()
            plt.ylabel(ylabs[i-1])
            i += 1
        for ts in sets:
            plt.subplot(2,3,i)
            plt.plot(t1.f, ts[0].data-ts[1].data, label=f'{labs[0]}-{labs[1]}')
            if len(ts)>2:
                plt.plot(t1.f, ts[1].data-ts[2].data, label=f'{labs[1]}-{labs[2]}')
            plt.legend()
            sl = slice(pGB["Frequency"]-1.5e-6, pGB["Frequency"]+1.5e-6)
            mima = np.max(np.abs(ts[1].sel(f=sl)-ts[2].sel(f=sl)))
            plt.axis([pGB["Frequency"]-1.5e-6, pGB["Frequency"]+1.5e-6, -mima, +mima])
            plt.grid()
            i += 1
        return tdi_fg
        
class SpritzViewer(PipeViewer):
    """
    """
    ymax = dict({"mbhb1":(4.1e-19,4.1e-19,4.1e-19,4.1e-19), 
                 "mbhb2":(4.1e-19,4e-21,4.1e-19,4e-21), 
                 "vgb":(2e-18,2e-18,4e-19,2e-22)})

    
    def __init__(self, path):
        self.path = path
        self.tdi2 = True
        
    def get_orbits(self):
        ofile = os.path.join(self.path, 'orbits.h5')
        orbits = Orbits.type(dict({'orbit_type':'file', 'filename':ofile, 'nominal_arm_length':2.5e9}))
        return orbits

    def load_tdi(self, source, noise=True, artifact=True):
        fn = f"{self.path}/{source}-training.h5"
        if noise and artifact:
            tdi = TDI.load(fn, name='obs/tdi')
        elif not noise and artifact:
            tdi = TDI.load(fn, name='noisefree/tdi')
        elif not artifact and noise:
            tdi = TDI.load(fn, name='clean/tdi')
        else:
            tdi = TDI.load(fn, name='sky/tdi')
        return tdi

    def get_all_dataset(self, source):
        mbhb1 = self.load_tdi(source, noise=True)
        mbhb1_nf = self.load_tdi(source, noise=False)
        mbhb1_af = self.load_tdi(source, artifact=False)
        mbhb1_gw = self.load_tdi(source, artifact=False, noise=False)
        ds = dict({"full":mbhb1, "noisefree":mbhb1_nf, "clean":mbhb1_af, "gwonly":mbhb1_gw})
        
        artifact_only = ds['full'] - ds["clean"]
        noise_only = ds["clean"] - ds["gwonly"]
        noise_free = ds["gwonly"] + artifact_only # alternative to mbhb_nf 

        
        ds["noiseonly"] = noise_only
        ds["artifactonly"] = artifact_only
        ds["noisefree2"] = noise_free
        return ds 
    
    def plot_all_dataset(self, source):

        ds = self.get_all_dataset(source)
        ymax = SpritzViewer.ymax
        
        plt.subplot(121)
        ds["full"]["X"].plot(label=f"{source} with noise and artifact")
        plt.legend(loc='lower left')
        plt.axis([None, None, -ymax[source][0], ymax[source][0]])
        plt.subplot(122)
        ds["noisefree"]["X"].plot(label=f"{source} without noise")
        plt.legend(loc='lower left')
        plt.axis([None, None, -ymax[source][1], ymax[source][1]])
        plt.figure(figsize=(15,5))
        plt.subplot(121)
        ds["clean"]["X"].plot(label=f"{source} without artifact")
        plt.legend(loc='lower left')
        plt.axis([None, None, -ymax[source][2], ymax[source][2]])
        plt.subplot(122)
        ds["gwonly"]["X"].plot(label=f"{source} without artifact and noise")
        plt.legend(loc='lower left')
        plt.axis([None, None, -ymax[source][3], ymax[source][3]])
        return ds


    def plot_combined_dataset(self, source):

        ds = self.get_all_dataset(source)
        ymax = SpritzViewer.ymax
        
        plt.figure(figsize=(15,5))
        plt.subplot(121)
        plt.plot(ds["full"].t, ds["noiseonly"]["X"], label=f"{source} clean - {source} gw only = noise only")
        plt.legend(loc='lower left')
        plt.axis([None, None, -ymax[source][2], ymax[source][2]])
        plt.subplot(122)
        plt.plot(ds["full"].t, ds["artifactonly"]["X"], label=f"{source} - {source} clean = artifact only")
        plt.legend(loc='lower left')
        plt.axis([None, None, -ymax[source][1], ymax[source][1]])
        return ds

    def compute_gwonly(self, source, from_y=True):

        ## Compute GW only TDI
        if from_y:
            tdi = self.load_tdi(source)
            return self.get_tdi_x_from_y(source, tdi.t.data)

        
        cat, units = self.get_catalog(source)
        if 'mbhb' in source:
            GW = HpHc.type("mbhb1", "MBHB", "IMRPhenomD")
            GW.set_param(dict(zip(cat.dtype.names, np.atleast_1d(cat)[0])))
        else:
            GW = HpHc.type("vgb", "GB", "TD_fdot")
            GW.set_param(cat)
            
        tdi = self.load_tdi(source)
        projector = ProjectedStrain(self.get_orbits())
        tvec = tdi.t.data
        t_min = tvec[0]
        t_max = tvec[-1]
        dt = tvec[1] - tvec[0]
        yArm = projector.arm_response(t_min, t_max, dt, GW.split())
        gw = TimeSeries(projector.compute_tdi_x(tvec, tdi2=self.tdi2), t0=t_min, dt=dt)
        return gw

   
    def plot_noise_psd(self, source):

        ds = self.get_all_dataset(source)
        
        f0, psdX0 = get_psd(ds["noiseonly"], k='X') #noise only
        f1, psdX1 = get_psd(ds["full"], k='X') #noise+gw+artifacts

        noise_model = "MRDv1"
        freq = np.logspace(-5, 1, 1000)
        Nmodel = get_noise_model(noise_model, freq)
        Npsd = Nmodel.psd(tdi2=True, option="X", freq=freq)
        
        plt.loglog(f1, np.sqrt(psdX1), label=source)
        plt.loglog(f0, np.sqrt(psdX0), label=f"{source}:noise only")
        plt.loglog(freq, np.sqrt(Npsd), color='k', label='MRDv1')
        plt.axis([2e-5, 0.3, None, None])
        plt.legend()

    def get_config(self, source):
        return h5io.load_config(os.path.join(self.path, f"{source}-training.h5"), "sky/config")
                                
    def get_pipe_config(self):
        source = 'mbhb1'
        cfg_pipe = h5io.load_config(os.path.join(self.path, f"{source}-training.h5"), "obs/config")
        cfg_pipe = dict([(k,v) for k,v in cfg_pipe.items() if 'seed' not in k]) # remove seed info
        return cfg_pipe
