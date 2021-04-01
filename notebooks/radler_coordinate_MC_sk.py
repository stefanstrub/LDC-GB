#%%
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
import scipy
import numpy as np
import xarray as xr
from astropy import units as u
import pandas as pd
import time
from copy import deepcopy

import torch
import gpytorch
from sklearn.metrics import mean_squared_error
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel, RBF

from botorch.models.gpytorch import GPyTorchModel
from botorch.acquisition.monte_carlo import qExpectedImprovement, qUpperConfidenceBound
from botorch.optim import optimize_acqf

import optuna

import ldc.io.hdf5 as hdfio
from ldc.lisa.noise import get_noise_model
from ldc.lisa.orbits import Orbits
from ldc.lisa.projection import ProjectedStrain
from ldc.common.series import TimeSeries, FrequencySeries, window
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr
from ldc.waveform.waveform import HpHc

from LISAhdf5 import LISAhdf5, ParsUnits
import tdi

import corner

# customized settings
plot_parameter = {  # 'backend': 'ps',
    "font.family": "serif",
    "font.serif": "times",
    "font.size": 16,
    "axes.labelsize": "medium",
    "axes.titlesize": "medium",
    "legend.fontsize": "medium",
    "xtick.labelsize": "medium",
    "ytick.labelsize": "medium",
    "grid.color": "k",
    "grid.linestyle": ":",
    "grid.linewidth": 0.5,
    "savefig.dpi": 150,
}

# tell matplotlib about your param_plots
rcParams.update(plot_parameter)
# set nice figure sizes
fig_width_pt = 464.0  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (np.sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
ratio = golden_mean
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
fig_width = fig_width_pt * inches_per_pt  # width in inches
fig_height = fig_width * ratio  # height in inches
fig_size = [fig_width, fig_height]
fig_size_squared = [fig_width, fig_width]
rcParams.update({"figure.figsize": fig_size})

# use a GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
dtype = torch.float


def Window(tm, offs=1000.0):

    xl = offs
    ind_r = np.argwhere(tm[-1] - tm <= offs)[0][0]
    xr = tm[ind_r]
    # print (xr)
    kap = 0.005
    winl = 0.5 * (1.0 + np.tanh(kap * (tm - xl)))
    winr = 0.5 * (1.0 - np.tanh(kap * (tm - xr)))
    # plt.plot(tm, winl)
    # plt.plot(tm, winr)
    # plt.grid(True)
    # plt.show()
    return winl * winr

np.random.seed(40)

parameters = [
    "Amplitude",
    "EclipticLatitude",
    "EclipticLongitude",
    "Frequency",
    "FrequencyDerivative",
    "Inclination",
    "InitialPhase",
    "Polarization",
]

DATAPATH = "/home/stefan/LDC/Radler/data"
sangria_fn = DATAPATH + "/dgb-tdi.h5"
sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
FD5 = LISAhdf5(sangria_fn)
Nsrc = FD5.getSourcesNum()
GWs = FD5.getSourcesName()
print("Found %d GW sources: " % Nsrc, GWs)
### TODO make sure GalBin is there
if GWs[0] != "GalBinaries":
    raise NotImplementedError
p = FD5.getSourceParameters(GWs[0])
td = FD5.getPreProcessTDI()
del_t = float(p.get("Cadence"))
Tobs = float(p.get("ObservationDuration"))

dt = del_t

# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = xr.Dataset(dict([(k, TimeSeries(td[:, n], dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))
# tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

# tdi_ts_mbhb = xr.Dataset(dict([(k,TimeSeries(tdi_ts_mbhb[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
# tdi_fs_mbhb = xr.Dataset(dict([(k,tdi_ts_mbhb[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_ts_dgb = xr.Dataset(dict([(k,TimeSeries(tdi_ts_dgb[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
# tdi_fs_dgb = xr.Dataset(dict([(k,tdi_ts_dgb[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_ts_igb = xr.Dataset(dict([(k,TimeSeries(tdi_ts_igb[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
# tdi_fs_igb = xr.Dataset(dict([(k,tdi_ts_igb[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_ts_vgb = xr.Dataset(dict([(k,TimeSeries(tdi_ts_vgb[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
# tdi_fs_vgb = xr.Dataset(dict([(k,tdi_ts_vgb[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
Npsd = Nmodel.psd()

# plt.figure(figsize=np.asarray(fig_size)*1.3)
# f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/10)
# # plt.loglog(tdi_fs["X"].f, (np.abs(tdi_fs["X"].values)**2/(len(tdi_fs["X"])*dt)), label="Data")
# plt.loglog(f, psdX, label="Data")
# f, psdX =  scipy.signal.welch(tdi_ts_vgb["X"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/10)
# plt.loglog(f, psdX, label="VGBs",zorder=5)
# f, psdX =  scipy.signal.welch(tdi_ts_mbhb["X"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/10)
# plt.loglog(f, psdX, label="MBHBs",zorder=4)
# f, psdX =  scipy.signal.welch(tdi_ts_dgb["X"].values+tdi_ts_igb["X"].values, fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/10)
# plt.loglog(f, psdX, label="GBs")
# plt.loglog(Nmodel.freq, Npsd, alpha=2, color='black', label='Instr. noise')
# plt.legend()
# plt.xlabel("Freq (Hz)")
# plt.ylabel("X-TDI PSD (1/Hz)")
# plt.grid(True)
# plt.axis([1e-5, 1e-1, 1e-43, 1e-34])
# plt.tight_layout()
# plt.show()

pGB = {}
ind = 0
for parameter in parameters:
    pGB[parameter] = p.get(parameter)[ind]

GB = fastGB.FastGB(delta_t=dt, T=float(tdi_ts["X"].t[-1]))  # in seconds

Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator="synthlisa")
fmin, fmax = float(Xs.f[0]), float(Xs.f[-1] + Xs.attrs["df"])
source = dict({"X": Xs, "Y": Ys, "Z": Zs})
start = time.time()
Xs_td, Ys_td, Zs_td = GB.get_td_tdixyz(template=pGB, simulator="synthlisa")


f_noise = np.logspace(-5, -1, 100)
Nmodel = get_noise_model(noise_model, f_noise)
freq = np.array(source["X"].sel(f=slice(fmin, fmax)).f)

Sn = Nmodel.psd(freq=freq, option="X")

######################################################
#%%
# find GB parameters
def loglikelihood(pGBs):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa")
    index_low = np.searchsorted(Xs.f, dataX.f[0])
    Xs = Xs[index_low : index_low + len(dataX)]
    Ys = Ys[index_low : index_low + len(dataY)]
    Zs = Zs[index_low : index_low + len(dataZ)]
    diff = np.abs(dataX - Xs.values) ** 2 + np.abs(dataY - Ys.values) ** 2 + np.abs(dataZ - Zs.values) ** 2
    # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
    p1 = float(np.sum(diff / Sn) * Xs.df) / 2.0
    # p1 = np.exp(p1)
    return -p1


def loglikelihoodratio(pGBs):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa")
    index_low = np.searchsorted(Xs.f, dataX.f[0])
    Xs = Xs[index_low : index_low + len(dataX)]
    Ys = Ys[index_low : index_low + len(dataY)]
    Zs = Zs[index_low : index_low + len(dataZ)]
    Xs.values = Xs - dataX.values
    Ys.values = Ys - dataY.values
    Zs.values = Zs - dataZ.values

    source = dict({"X": Xs, "Y": Ys, "Z": Zs})
    data = dict({"X": dataX, "Y": dataY, "Z": dataZ})
    # SNRd = compute_tdi_snr(source, Nmodel, data=data)["tot2"]
    SNRs = compute_tdi_snr(source, Nmodel)["tot2"]
    # p1 = SNRd - 0.5*SNRs
    # diff = np.abs(dataX - Xs.values)**2 + np.abs(dataY - Ys.values)**2 + np.abs(dataZ - Zs.values)**2
    # # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
    # p1 = float(np.sum(diff / (Sn+noise))*Xs.df)/2.0
    # p1 = np.exp(p1)
    return SNRs


# Number of histogram bins.
n_bin = 50
# Total number of proposed samples.
number_of_samples = 8 * 10 ** 1
cutoff_ratio = 1000

# Neil2005
# pGB['Amplitude'] = 1.4e-22
# pGB['EclipticLatitude'] = 0.399
# pGB['EclipticLongitude'] = 5.71
# pGB['Frequency'] = 1.0020802e-3
# pGB['FrequencyDerivative'] = 1e-19
# pGB['Inclination'] = 0.96
# pGB['InitialPhase'] = 1.0
# pGB['Polarization'] = 1.3

parametersfd = [
    "Amplitude",
    "EclipticLatitude",
    "EclipticLongitude",
    "Frequency",
    "FrequencyDerivative",
    "Inclination",
    "InitialPhase",
    "Polarization",
]
boundaries = {
    "Amplitude": [10 ** -23.0, 2 * 10 ** -22.0],
    "EclipticLatitude": [-1.0, 1.0],
    "EclipticLongitude": [-np.pi, np.pi],
    "Frequency": [pGB["Frequency"] * 0.9999, pGB["Frequency"] * 1.0001],
    "FrequencyDerivative": [10 ** -20.0, 10 ** -16.0],
    "Inclination": [-1.0, 1.0],
    "InitialPhase": [0.0, 1.0 * np.pi],
    "Polarization": [np.pi, 2.0 * np.pi],
}

previous_max = [0.2090, 0.1000, 0.8469, 0.5276, 0.7168, 0.9667, 0.0970, 0.0000]
# previous_max = [0.2090, 0.2667, 0.7333, 0.5276, 0.7168, 0.9667, 0.0970, 0.0000]
# previous_max = [0.2090, 0.3333, 0.7333, 0.5276, 0.7168, 0.8667, 0.0970, 0.0000]
# previous_max = [0.2090, 0.3667, 0.6667, 0.5276, 0.7168, 0.8667, 0.0970, 0.0000]
# previous_max = [0.4000, 0.3667, 0.6667, 0.5276, 0.7168, 0.8667, 0.0970, 0.0000]
# previous_max = [0.4000, 0.3667, 0.6667, 0.5276, 0.7667, 0.8667, 0.0970, 0.0000]
# previous_max = [0.2333, 0.3667, 0.6667, 0.5276, 0.9667, 0.8667, 0.0970, 0.0000]
# previous_max = [0.333,0.54,0.134,0.7345,0.6456,0.2645,0.8216,0.000]
# previous_max = [0.2090, 0.1000, 0.8469, 0.3276, 0.7168, 0.9667, 0.0970, 0.0000]
# previous_max = [0.45090, 0.5600, 0.123469, 0.87276, 0.2341168, 0.56667, 0.5689970, 0.0000]
# previous_max = [0.45090, 0.5600, 0.123469, 0.3276, 0.2341168, 0.56667, 0.9689970, 0.0000]
# previous_max = [0.86436456, 0.3156825,  0.6350386,  0.55715334, 0.5604474,  0.7789818, 0.03608589, 0.0]
# previous_max =[0.80888367, 0.35581076, 0.62365836, 0.5551591,  0.5607991,  0.76172084,0.03608589, 0.0]
# previous_max =[0.80888367, 0.35581076, 0.62365836, 0.5,  0.5607991,  0.76172084,0.03608589, 0.0]
previous_max = np.random.rand(8)
i = 0
pGBs = deepcopy(pGB)
for parameter in parameters:
    if parameter in ["FrequencyDerivative"]:
        i -= 1
    elif parameter in ["EclipticLatitude"]:
        pGBs[parameter] = np.arcsin((previous_max[i] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
    elif parameter in ["Inclination"]:
        pGBs[parameter] = np.arccos((previous_max[i] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
    else:
        pGBs[parameter] = (previous_max[i] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
        print(parameter, pGBs[parameter], previous_max[i])
    i += 1


Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa")
# dataX2,dataY2,dataZ2 = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator='synthlisa')
psd_signal = np.abs(Xs.values) ** 2 + np.abs(Ys.values) ** 2 + np.abs(Zs.values) ** 2
highSNR = psd_signal > np.max(psd_signal) / cutoff_ratio
lowerindex = np.where(highSNR)[0][0] - 10
higherindex = np.where(highSNR)[0][-1] + 10
dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin + len(Xs)))[lowerindex:higherindex]
dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin + len(Ys)))[lowerindex:higherindex]
dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin + len(Zs)))[lowerindex:higherindex]
Xs, Ys, Zs = (
    Xs[lowerindex:higherindex],
    Ys[lowerindex:higherindex],
    Zs[lowerindex:higherindex],
)
# index_low = np.searchsorted(dataX2.f, Xs.f[0])
# dataX = dataX2[index_low:index_low+len(Xs)]
# dataY = dataY2[index_low:index_low+len(Xs)]
# dataZ = dataZ2[index_low:index_low+len(Xs)]
spd_data = np.abs(dataX) ** 2 + np.abs(dataY) ** 2 + np.abs(dataZ) ** 2
noise = (np.mean(spd_data[:2]) + np.mean(spd_data[-2:])).values / 2
noise = 0  # (np.mean(spd_data).values)/2
fmin, fmax = float(Xs.f[0]), float(Xs.f[-1] + Xs.attrs["df"])
freq = np.array(Xs.sel(f=slice(fmin, fmax)).f)
Nmodel = get_noise_model(noise_model, freq)
Sn = Nmodel.psd(freq=freq, option="X")
diffd = np.abs(dataX) ** 2 + np.abs(dataY) ** 2 + np.abs(dataZ) ** 2
pd = float(np.sum(diffd / (Sn + noise)) * Xs.df) / 2.0
diff = np.abs(dataX - Xs.values) ** 2 + np.abs(dataY - Ys.values) ** 2 + np.abs(dataZ - Zs.values) ** 2
p1 = float(np.sum(diff / (Sn + noise)) * Xs.df) / 2.0
p1 = -p1
# p1 = np.exp(p1)

# def loglikelihoodstas(pGBs):
#     Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator='synthlisa')
#     index_low = np.searchsorted(Xs.f, dataX.f[0])
#     Xs = Xs[index_low:index_low+len(dataX)]
#     Ys = Ys[index_low:index_low+len(dataY)]
#     Zs = Zs[index_low:index_low+len(dataZ)]
#     # fr = np.arange(Af.kmin, Af.kmin+len(Af))*df
#     ### TODO I assume that the frequency range is the same for A and E templates

# #     SA = tdi.noisepsd_AE(fr, model='Proposal', includewd=None)
#     SA = tdi.noisepsd(Xs.f, model='Proposal', includewd=None)
#     SNR2 = np.sum( np.real(DAf[ib:ie] * np.conjugate(Af.data) + DEf[ib:ie] * np.conjugate(Ef.data))/SA )
#     hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /SA)

#     print ("SN", 4.0*df*SNR2, 4.0*df* hh)

#     loglik = 4.0*df*( SNR2 - 0.5 * hh )
#     return (loglik)

# start = time.time()
# for i in range(10**2):
#     p = loglikelihood(pGB)
# print(time.time()-start, p)
# start = time.time()
# for i in range(10**2):
#     p = loglikelihoodratio(pGB)
# print(time.time()-start, p)
# for i in range(10**2):
#     p = loglikelihoodstas(pGB)
# print(time.time()-start, p)

# frequency_lower_boundary = dataX.f[0].values+(dataX.f[-1].values - dataX.f[0].values)*4/10
# frequency_upper_boundary = dataX.f[-1].values-(dataX.f[-1].values - dataX.f[0].values)*3/10
# [0.0004725, 0.0004727]
# boundaries_small = deepcopy(boundaries)
# part_ratio = 10
# for parameter in parameters:
#     if parameter in ['EclipticLongitude','Frequency']:
#         boundaries_small[parameter] = [pGB[parameter]-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,pGB[parameter]+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
#     if parameter == 'EclipticLatitude':
#         boundaries_small[parameter] = [np.sin(pGB[parameter])-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,np.sin(pGB[parameter])+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
#     elif parameter == 'Inclination':
#         boundaries_small[parameter] = [np.cos(pGB[parameter])-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,np.cos(pGB[parameter])+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
#     else:
#         boundaries_small[parameter] = [pGB[parameter]-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,pGB[parameter]+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
# boundaries = boundaries_small

Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator="synthlisa")
index_low = np.searchsorted(Xs.f, dataX.f[0])
Xs = Xs[index_low : index_low + len(dataX)]
Ys = Ys[index_low : index_low + len(dataY)]
Zs = Zs[index_low : index_low + len(dataZ)]

# plt.figure(figsize=fig_size)
# ax1 = plt.subplot(111)
# # plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
# ax1.plot(dataX.f * 1000, dataX.values, label="data", marker=".", zorder=5)
# ax1.plot(Xs.f * 1000, Xs.values, label="VGB", marker=".", zorder=5)
# ax1.plot(
#     Xs.f * 1000,
#     dataX.values - Xs.values,
#     label="residual",
#     alpha=0.8,
#     color="red",
#     marker=".",
# )
# plt.legend()
# plt.show()

# plt.figure()
# ax1=plt.subplot(111)
# # plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
# ax1.loglog(tdi_fs['X'].f*1000,np.abs(tdi_fs['X'].values), label='data',zorder=1)
# ax1.loglog(Xs.f*1000, np.abs(Xs.values), label='VGB',zorder=5)
# # ax1.loglog(Xs.f*1000, np.abs(dataX.values-Xs.values), label='residual',alpha=0.8, color='red',marker='.')
# ax1.loglog(tdi_fs_mbhb['X'].f*1000, np.abs(tdi_fs_mbhb['X']), label='MBHBs',alpha=0.8,color='green')
# ax1.loglog(tdi_fs_dgb['X'].f*1000, np.abs(tdi_fs_dgb['X'].values+tdi_fs_igb['X'].values), label='GBs',alpha=0.8,color='red')

# ax1.set_xlabel('f [mHz]')
# ax1.set_ylabel('X-TDI [1/Hz]')
# ax1.legend()
# plt.show()

samples = xr.Dataset(
    dict(
        [
            (
                name,
                xr.DataArray(
                    np.zeros(number_of_samples),
                    dims=("number_of_sample"),
                    coords={"number_of_sample": range(number_of_samples)},
                ),
            )
            for name, titles in pGBs.items()
        ]
    )
)
pGBs01 = {}
for parameter in parameters:
    if parameter in ["EclipticLatitude"]:
        samples[parameter][0] = (np.sin(pGBs[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
    elif parameter in ["Inclination"]:
        samples[parameter][0] = (np.cos(pGBs[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
    else:
        samples[parameter][0] = (pGBs[parameter] - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
    pGBs01[parameter] = samples[parameter][0]
samples["Likelihood"] = samples["Inclination"]
samples["Likelihood"][0] = p1

print("p start", p1, "p true", loglikelihood(pGB))


def sampler(
    number_of_samples,
    parameters,
    pGB,
    boundaries,
    p1,
    uniform=False,
    MCMC=False,
    only=False,
    onlyparameter="Frequency",
    twoD=False,
    secondparameter="Amplitude",
    calculate_loglikelihood=True,
):
    samples = xr.Dataset(
        dict(
            [
                (
                    name,
                    xr.DataArray(
                        np.zeros(number_of_samples),
                        dims=("number_of_sample"),
                        coords={"number_of_sample": range(number_of_samples)},
                    ),
                )
                for name, titles in pGB.items()
            ]
        )
    )
    samples = {}
    pGBs01 = {}
    for parameter in parameters:
        samples[parameter] = []
        if parameter in ["EclipticLatitude"]:
            samples[parameter].append((np.sin(pGB[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0]))
        elif parameter in ["Inclination"]:
            samples[parameter].append((np.cos(pGB[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0]))
        else:
            samples[parameter].append((pGB[parameter] - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0]))
        pGBs01[parameter] = samples[parameter][0]
    samples["Likelihood"] = []
    psample = loglikelihood(pGB)
    samples["Likelihood"].append(psample)
    j = 0
    number_of_sampels_sqrt = np.sqrt(number_of_samples)
    for i in range(1, number_of_samples):
        if only:
            parameter = onlyparameter
            if uniform:
                pGBs01[parameter] = i / number_of_samples
            else:
                pGBs01[parameter] = np.random.rand()
            if parameter != "InitialPhase":
                pGBs01["InitialPhase"] = samples["InitialPhase"][0]
            if parameter != "Polarization":
                pGBs01["Polarization"] = samples["Polarization"][0]
        elif twoD:
            parameter = onlyparameter
            parameter2 = secondparameter
            if uniform:
                if i % number_of_sampels_sqrt == 0:
                    j += 1
                pGBs01[parameter] = ((i - 1) % number_of_sampels_sqrt + 1) / (number_of_sampels_sqrt + 2)
                pGBs01[parameter2] = (j + 1) / (number_of_sampels_sqrt + 2)
            else:
                pGBs01[parameter] = np.random.rand()
                pGBs01[parameter2] = np.random.rand()
            # if parameter != 'InitialPhase':
            #     pGBs01['InitialPhase'] = samples['InitialPhase'][0]
        else:
            for parameter in parameters:
                pGBs01[parameter] = np.random.rand()
                # if parameter in ['Amplitude']:#,'FrequencyDerivative','Amplitude','EclipticLongitude']:
                #     pGBs01[parameter] = i/number_of_samples
                # elif parameter in ['FrequencyDerivative']:
                #     pass
                # if parameter in ['Inclination']:
                #     pGBs01[parameter] = np.random.rand()
                # elif parameter in ['Inclination']:
                #     pGBs01[parameter] = np.random.rand()
                # else:
                #     pGBs01[parameter] = np.random.rand()
        for parameter in parameters:
            # if parameter in ['FrequencyDerivative']:
            #     pass
            if parameter in ["EclipticLatitude"]:
                pGBs[parameter] = np.arcsin((pGBs01[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
            elif parameter in ["Inclination"]:
                pGBs[parameter] = np.arccos((pGBs01[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
            else:
                pGBs[parameter] = (pGBs01[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
        # pGBs01array = np.zeros(len(parametersfd))
        # i = 0
        # for name in parametersfd:
        #     pGBs01array[i] = pGBs01[name]
        #     i +=1
        # pGBscaled = scaletooriginal(pGBs01array ,boundaries)
        # print(loglikelihood(pGBs), loglikelihood(pGBscaled))
        if calculate_loglikelihood:
            psample = loglikelihood(pGBs)
            samples["Likelihood"].append(psample)
        for parameter in parameters:
            samples[parameter].append(pGBs01[parameter])
    for parameter in parameters:
        samples[parameter] = np.asarray(samples[parameter])
    samples["Likelihood"] = np.asarray(samples["Likelihood"])
    return samples


#%%


class ExactGPModel(gpytorch.models.ExactGP, GPyTorchModel):
    _num_outputs = 1  # to inform GPyTorchModel API

    def __init__(self, train_x, train_y, likelihood, kernel):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ZeroMean()
        kernelA = gpytorch.kernels.RBFKernel(
            active_dims=(0),
            lengthscale_constraint=gpytorch.constraints.GreaterThan(0.2),
        )
        kernelLat = gpytorch.kernels.PeriodicKernel(
            active_dims=torch.tensor([1]),
            lengthscale_constraint=gpytorch.constraints.GreaterThan(0.8),
        )
        kernelLong = gpytorch.kernels.PeriodicKernel(
            active_dims=torch.tensor([2]),
            lengthscale_constraint=gpytorch.constraints.GreaterThan(0.8),
        )
        kernelF = gpytorch.kernels.RBFKernel(
            active_dims=(3),
            lengthscale_constraint=gpytorch.constraints.Interval(0.06, 0.1),
        )
        # kernelFD = gpytorch.kernels.RBFKernel(active_dims=(4))
        kernelI = gpytorch.kernels.PolynomialKernel(power=2, active_dims=(4))
        # kernelIP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([6]))
        kernelIP = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.CosineKernel(
                active_dims=torch.tensor([5]),
                period_length_constraint=gpytorch.constraints.Interval(0.499, 0.501),
            ),
            outputscale_constraint=gpytorch.constraints.Interval(0.9, 1.1),
        )
        # kernelP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([7]))
        # kernelP = gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]))*gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]))
        kernelP = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.CosineKernel(
                active_dims=torch.tensor([6]),
                period_length_constraint=gpytorch.constraints.Interval(0.249, 0.251),
            ),
            outputscale_constraint=gpytorch.constraints.Interval(0.9, 1.1),
        )
        # kernel = kernelA + kernelLat + kernelLong + kernelF + kernelFD + kernelI + kernelIP + kernelP
        kernel = kernelA * kernelLat * kernelLong * kernelF * kernelI * kernelIP * kernelP
        # self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel(ard_num_dims=8))
        self.covar_module = gpytorch.kernels.RBFKernel(ard_num_dims=8)
        self.to(train_x)  # make sure we're on the right device/dtype

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


number_of_test_samples = 20
test_x = {}
test_y = {}


def objective(trial):
    for parameter in changeableparameters:
        parametervalue = trial.suggest_uniform(parameter, boundaries[parameter][0], boundaries[parameter][1])
        maxpGB2[parameter] = parametervalue
        if parameter in ["EclipticLatitude"]:
            maxpGB2[parameter] = np.arcsin(parametervalue)
        elif parameter in ["Inclination"]:
            maxpGB2[parameter] = np.arccos(parametervalue)
    p = loglikelihood(maxpGB2)
    return p


def model(optimize_parameters):
    i = 0
    for parameter in changeableparameters:
        parametervalue = optimize_parameters[i]
        maxpGB2[parameter] = (parametervalue * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
        if parameter in ["EclipticLatitude"]:
            maxpGB2[parameter] = troch.arcsin(parametervalue)
        elif parameter in ["Inclination"]:
            maxpGB2[parameter] = torch.arccos(parametervalue)
    p = loglikelihood(maxpGB2)
    i += 1
    return p


def scaletooriginal(previous_max, boundaries):
    maxpGB = deepcopy(pGBs)
    for parameter in parametersfd:
        if parameter in ["EclipticLatitude"]:
            maxpGB[parameter] = np.arcsin((previous_max[parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        elif parameter in ["Inclination"]:
            maxpGB[parameter] = np.arccos((previous_max[parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        else:
            maxpGB[parameter] = (previous_max[parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
    return maxpGB


def plotplanes(parameterstocheck, parameter2, plot_x, plot_y):
    fig, ax = plt.subplots(2, 4, figsize=(15, 15))
    plt.suptitle("loglikelihood true")
    i = 0
    for parameter in parameterstocheck:
        if parameter != parameter2:
            j = 0
            if i > 3:
                j = 1
            with torch.no_grad():
                ax[j, i % 4].axvline(x=previous_max[parametersfd.index(parameter)], color="k")
                ax[j, i % 4].axhline(y=previous_max[parametersfd.index(parameter2)], color="k")
                im = ax[j, i % 4].scatter(
                    plot_x[parameter + parameter2].numpy()[:, parametersfd.index(parameter)],
                    plot_x[parameter + parameter2].numpy()[:, parametersfd.index(parameter2)],
                    c=plot_y[parameter + parameter2].numpy()[:],
                )
                ax[j, i % 4].set_xlabel(parameter)
                ax[j, i % 4].set_ylabel(parameter2)
                fig.colorbar(im, ax=ax[j, i % 4])
        i += 1
    plt.show()

def traingpmodelsk(train_x, train_y, kernel, sigma, nu):
    train_x = train_x.numpy()
    train_y = train_y.numpy()
    kernel = RBF(length_scale=[1,1,1,1,1,1,1,1])
    gpr = GaussianProcessRegressor(kernel=kernel,
            random_state=0).fit(train_x, train_y)
    gpr.score(train_x, train_y)
    return gpr
def traingpmodel(train_x, train_y, kernel, sigma, nu):
    bo_iterations = 1
    for bo_iter in range(bo_iterations):
        # initialize likelihood and model
        likelihood = gpytorch.likelihoods.GaussianLikelihood()
        model = ExactGPModel(train_x, train_y, likelihood, kernel)

        training_iter = 15

        hypers = {
            "likelihood.noise": torch.tensor(0.0001),
        }
        model.initialize(**hypers)
        # Find optimal model hyperparameters
        model.train()
        likelihood.train()

        # Use the adam optimizer
        optimizer = torch.optim.Adam(
            [
                # {"params": model.mean_module.parameters()},
                {"params": model.covar_module.parameters()},
            ],
            lr=0.1,
        )  # Includes GaussianLikelihood parameters
        # optimizer = torch.optim.Adam(model.parameters(), lr=0.1)

        # "Loss" for GPs - the marginal log likelihood
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

        for i in range(training_iter):
            torch.cuda.empty_cache()
            # Zero gradients from previous iteration
            optimizer.zero_grad()
            # Output from model
            output = model(train_x)
            # Calc loss and backprop gradients
            loss = -mll(output, train_y)
            loss.backward()
            print(
                "Iter %d/%d - Loss: %.3f  Len: %.3f Len: %.3f Len: %.3f Len: %.3f Len: %.3f Len: %.3f Len: %.3f Len: %.3f"
                % (
                    i + 1,
                    training_iter,
                    loss.item(),
                    model.covar_module.lengthscale[0][0],
                    model.covar_module.lengthscale[0][1],
                    model.covar_module.lengthscale[0][2],
                    model.covar_module.lengthscale[0][3],
                    model.covar_module.lengthscale[0][4],
                    model.covar_module.lengthscale[0][5],
                    model.covar_module.lengthscale[0][6],
                    model.covar_module.lengthscale[0][7],
                )
            )
            optimizer.step()

        best_value = train_y.max()
        print("best value", bo_iter, (best_value * sigma) + nu, train_x[train_y.argmax()])
        best_value = best_value
        if bo_iter < bo_iterations - 1:
            qEI = qExpectedImprovement(model=model, best_f=best_value)
            qUCB = qUpperConfidenceBound(model=model, beta=0.001)
            UCB = UpperConfidenceBound(model=model, beta=0.001, maximize=True)
            candidates = 1
            new_point_analytic, _ = optimize_acqf(
                acq_function=UCB,
                bounds=torch.tensor([[0.0] * 8, [1.0] * 8]),
                q=candidates,
                num_restarts=1,
                raw_samples=100,
                options={},
            )
            train_x = torch.cat((train_x, new_point_analytic), 0)
            new_point_analytic = new_point_analytic.cpu().numpy()
            for candi_n in range(candidates):
                para_n = 0
                for parameter in parameters:
                    if parameter in ["EclipticLatitude"]:
                        pGBs[parameter] = np.arcsin((new_point_analytic[candi_n][para_n] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
                    elif parameter in ["Inclination"]:
                        pGBs[parameter] = np.arccos((new_point_analytic[candi_n][para_n] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
                    else:
                        pGBs[parameter] = (new_point_analytic[candi_n][para_n] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
                    para_n += 1
                loglike = loglikelihood(pGBs)
                print("new point", i, new_point_analytic[candi_n], loglike)

                loglike = (loglike - nu) / sigma
                loglike = torch.tensor([loglike]).float()
                train_y = torch.cat((train_y, loglike), 0)
    return model, likelihood


maxpGB = deepcopy(pGBs)
notreduced = False
notreduced2 = False
boundaries_reduced = deepcopy(boundaries)
best_params = deepcopy(pGBs)
best_value = p1
parameters_recorded = []
best_run = 0
for n in range(15):
    parameters_recorded1 = {}
    no_improvement_counter = 0
    for parameter in parametersfd:
        parametervalue = np.random.uniform(boundaries[parameter][0], boundaries[parameter][1])
        maxpGB[parameter] = parametervalue
        if parameter in ["EclipticLatitude"]:
            maxpGB[parameter] = np.arcsin(parametervalue)
        elif parameter in ["Inclination"]:
            maxpGB[parameter] = np.arccos(parametervalue)
        parameters_recorded1[parameter] = []
        parameters_recorded1[parameter].append(maxpGB[parameter])
    parameters_recorded1["Loglikelihood"] = []
    if n == 0:
        maxpGB = deepcopy(pGBs)
    previous_best = loglikelihood(maxpGB)
    maxpGB2 = deepcopy(maxpGB)
    for i in range(100):
        parameter1 = parametersfd[i % 7]
        parameter2 = parametersfd[np.random.randint(0, 6)]
        parameter3 = parametersfd[np.random.randint(0, 6)]
        # parameter2 = 'InitialPhase'
        # parameter1 = 'Inclination'
        while parameter2 == parameter1:
            parameter2 = parametersfd[np.random.randint(0, 6)]
        while parameter3 == parameter1 or parameter3 == parameter2:
            parameter3 = parametersfd[np.random.randint(0, 6)]
        # if parameter1 == 'Frequency':
        #     parameter2 = 'Polarization'
        parametersreduced = [parameter1]
        changeableparameters = [parameter1, parameter2, parameter3]
        params = np.zeros(len(changeableparameters))

        optuna.logging.set_verbosity(optuna.logging.WARNING)
        start = time.time()
        study = optuna.create_study(sampler=optuna.samplers.RandomSampler(), direction="maximize")
        study.optimize(objective, n_trials=50)
        print("optuna time", time.time() - start)
        if study.best_value > previous_best:
            no_improvement_counter = 0
            previous_best = study.best_value
            for parameter in changeableparameters:
                maxpGB[parameter] = study.best_params[parameter]
                if parameter in ["EclipticLatitude"]:
                    maxpGB[parameter] = np.arcsin(study.best_params[parameter])
                elif parameter in ["Inclination"]:
                    maxpGB[parameter] = np.arccos(study.best_params[parameter])
        else:
            no_improvement_counter += 1
        if no_improvement_counter > 16:
            past_mean = 0
            sum_count = 0
            for l in range(n):
                s = 0
                try:
                    past_mean += parameters_recorded[l]["Loglikelihood"][i]
                    sum_count += 1
                except:
                    pass
            try:
                past_mean = past_mean / sum_count
                if previous_best > past_mean:
                    no_improvement_counter = 0
                else:
                    print("no improvement")
                    break
            except:
                no_improvement_counter = 0
                pass
             
        if i in [20,30,40,50]:
            past_mean = 0
            sum_count = 0
            for l in range(n):
                s = 0
                try:
                    past_mean += parameters_recorded[l]["Loglikelihood"][i]
                    sum_count += 1
                except:
                    pass
            try:
                past_mean = past_mean / sum_count
                if previous_best > past_mean:
                    pass
                else:
                    break
            except:
                pass

        start = time.time()
        print(i, previous_best, loglikelihood(maxpGB), maxpGB)
        parameters_recorded1["Loglikelihood"].append(loglikelihood(maxpGB))
        for parameter in parametersfd:
            parameters_recorded1[parameter].append(maxpGB[parameter])

        maxpGB2 = deepcopy(maxpGB)
        x = 1
    if previous_best > best_value:
        best_run = n
        best_value = previous_best
        best_params = deepcopy(maxpGB)
    parameters_recorded.append(parameters_recorded1)
    # previous_max, max_loglike, observed_pred = planeAdam(maxpGB,parametersreduced,parameter2,resolution, boundaries)
    # maxpGB = scaletooriginal(previous_max,boundaries)
    # print('random time',time.time()-start)
    # print(i,max_loglike, maxpGB)
    # print(previous_max)
    # real_loglike = loglikelihood(maxpGB)
    # plotplanes(parametersreduced,parameter2, test_x, test_y)
    # if parameter1 == 'Frequency':
    # plotplanes(parametersreduced,parameter2, test_x, observed_pred)
    # if real_loglike > -0.9 and notreduced:
    #     notreduced = False
    #     ratio = 0.2
    #     for parameter in parameters:
    #         length = boundaries[parameter][1]-boundaries[parameter][0]
    #         if parameter == 'EclipticLatitude':
    #             boundaries_reduced[parameter] = [np.sin(maxpGB[parameter])-length*ratio/2, np.sin(maxpGB[parameter])+length*ratio/2]
    #         elif parameter == 'Inclination':
    #             boundaries_reduced[parameter] = [np.cos(maxpGB[parameter])-length*ratio/2, np.cos(maxpGB[parameter])+length*ratio/2]
    #         else:
    #             boundaries_reduced[parameter] = [maxpGB[parameter]-length*ratio/2, maxpGB[parameter]+length*ratio/2]
    #         if boundaries_reduced[parameter][0] <  boundaries[parameter][0]:
    #             boundaries_reduced[parameter][0] =  boundaries[parameter][0]
    #         if boundaries_reduced[parameter][1] >  boundaries[parameter][1]:
    #             boundaries_reduced[parameter][1] =  boundaries[parameter][1]
    #     boundaries = boundaries_reduced
    # if real_loglike > -0.8 and notreduced2:
    #     notreduced2 = False
    #     ratio = 0.5
    #     for parameter in parameters:
    #         length = boundaries[parameter][1]-boundaries[parameter][0]
    #         if parameter == 'EclipticLatitude':
    #             boundaries_reduced[parameter] = [np.sin(maxpGB[parameter])-length*ratio/2, np.sin(maxpGB[parameter])+length*ratio/2]
    #         elif parameter == 'Inclination':
    #             boundaries_reduced[parameter] = [np.cos(maxpGB[parameter])-length*ratio/2, np.cos(maxpGB[parameter])+length*ratio/2]
    #         else:
    #             boundaries_reduced[parameter] = [maxpGB[parameter]-length*ratio/2, maxpGB[parameter]+length*ratio/2]
    #         if boundaries_reduced[parameter][0] <  boundaries[parameter][0]:
    #             boundaries_reduced[parameter][0] =  boundaries[parameter][0]
    #         if boundaries_reduced[parameter][1] >  boundaries[parameter][1]:
    #             boundaries_reduced[parameter][1] =  boundaries[parameter][1]
    #     boundaries = boundaries_reduced

maxpGB = best_params
ratio = 0.1
for parameter in parameters:
    length = boundaries[parameter][1] - boundaries[parameter][0]
    if parameter == "EclipticLatitude":
        boundaries_reduced[parameter] = [
            np.sin(maxpGB[parameter]) - length * ratio / 2,
            np.sin(maxpGB[parameter]) + length * ratio / 2,
        ]
    elif parameter == "Inclination":
        boundaries_reduced[parameter] = [
            np.cos(maxpGB[parameter]) - length * ratio / 2,
            np.cos(maxpGB[parameter]) + length * ratio / 2,
        ]
    elif parameter == "Frequency":
        boundaries_reduced[parameter] = [
            maxpGB[parameter] - length * ratio / 8,
            maxpGB[parameter] + length * ratio / 8,
        ]
    elif parameter == "Amplitude":
        boundaries_reduced[parameter] = [
            maxpGB[parameter] - length * ratio / 2,
            maxpGB[parameter] + length * ratio / 2,
        ]
    # elif parameter in ["InitialPhase",'Polarization']:
    #     boundaries_reduced[parameter] = [
    #         pGB[parameter] - length * ratio / 2*4,
    #         pGB[parameter] + length * ratio / 2*4,
    #     ]
    #     boundaries_reduced[parameter] = boundaries[parameter]
    else:
        boundaries_reduced[parameter] = [
            maxpGB[parameter] - length * ratio / 2,
            maxpGB[parameter] + length * ratio / 2,
        ]
    if boundaries_reduced[parameter][0] < boundaries[parameter][0]:
        boundaries_reduced[parameter][0] = boundaries[parameter][0]
    if boundaries_reduced[parameter][1] > boundaries[parameter][1]:
        boundaries_reduced[parameter][1] = boundaries[parameter][1]

resolution = 2000
# if parameter != parameter2:
resolution_reduced = int(20 ** 2)
resolution_reduced = resolution
# if 'Frequency' in [parameter,parameter2]:
#     resolution_reduced = int(15**2)
start = time.time()
parameter = "Frequency"
train_samples = sampler(
    resolution_reduced,
    parameters,
    maxpGB,
    boundaries_reduced,
    p1,
    uniform=False,
    twoD=False,
    onlyparameter=parameter,
    secondparameter=parameter2,
)
print(time.time() - start)
train_x = np.zeros((resolution_reduced, len(parameters)))
i = 0
for name in parametersfd:
    train_x[:, i] = train_samples[name]
    i += 1
train_y = train_samples["Likelihood"]
train_x = torch.from_numpy(train_x).float()
train_y = torch.from_numpy(train_y).float()
# if parameter != parameter2:
resolution = 1000
test_samples = sampler(
    resolution,
    parameters,
    maxpGB,
    boundaries_reduced,
    p1,
    uniform=False,
    twoD=False,
    onlyparameter=parameter,
    secondparameter=parameter2,
    calculate_loglikelihood=True,
)
test_x[parameter + parameter2] = np.zeros((resolution, len(parameters)))
i = 0
for name in parametersfd:
    test_x[parameter + parameter2][:, i] = test_samples[name]
    i += 1
test_y[parameter + parameter2] = test_samples["Likelihood"]
test_x[parameter + parameter2] = torch.from_numpy(test_x[parameter + parameter2]).float()
test_y[parameter + parameter2] = torch.from_numpy(test_y[parameter + parameter2]).float()

kernel = gpytorch.kernels.RBFKernel(ard_num_dims=2)

# train_y = np.exp(train_y-loglikelihood(maxpGB))
nu = np.mean(train_y.numpy())
sigma = np.std(train_y.numpy())
train_y = (train_y - nu) / sigma
model, likelihood = traingpmodel(train_x, train_y, kernel, sigma, nu)
start = time.time()
gpr = traingpmodelsk(train_x, train_y, kernel, sigma, nu)
print(time.time() - start)
model.eval()
likelihood.eval()
start = time.time()
observed_pred_sk = gpr.predict(test_x[parameter + parameter2])
print(time.time() - start)
observed_pred_sk_scaled = observed_pred_sk*sigma +nu
print("sqrt(MSE) ",parameter + parameter2,np.sqrt(mean_squared_error(test_y[parameter + parameter2].numpy(),observed_pred_sk_scaled)))
observed_pred = {}
observed_pred_mean = {}
with torch.no_grad(), gpytorch.settings.fast_pred_var():
    observed_pred[parameter + parameter2] = likelihood(model(test_x[parameter + parameter2]))
    observed_pred_mean[parameter + parameter2] = (observed_pred[parameter + parameter2].mean * sigma) + nu
print("sqrt(MSE) ",parameter + parameter2,np.sqrt(mean_squared_error(test_y[parameter + parameter2].numpy(),observed_pred_mean[parameter + parameter2].numpy())))

resolution = 10 ** 7
start = time.time()
test_samples = sampler(
    resolution,
    parameters,
    maxpGB,
    boundaries_reduced,
    p1,
    uniform=False,
    twoD=False,
    onlyparameter=parameter,
    secondparameter=parameter2,
    calculate_loglikelihood=False,
)
print('sample time', time.time()-start)
test_x_m = np.zeros((resolution, len(parameters)))
i = 0
for name in parametersfd:
    test_x_m[:, i] = test_samples[name]
    i += 1
test_x_m = torch.from_numpy(test_x_m).float()
start = time.time()
observed_pred_sk = gpr.predict(test_x_m[:5*10**5])
for i in range(int(resolution/(5*10**5))-1):
    observed_pred_sk = np.append(observed_pred_sk,gpr.predict(test_x_m[(i+1)*5*10**5:(i+2)*5*10**5]))
# observed_pred_sk = gpr.predict(test_x_m[:5*10**5])
# observed_pred_sk = gpr.predict(test_x_m)
observed_pred_mean = observed_pred_sk*sigma +nu
print('eval time', time.time()-start)
# with torch.no_grad(), gpytorch.settings.fast_pred_var():
#     observed_pred = likelihood(model(test_x_m))
#     observed_pred_mean = (observed_pred.mean * sigma) + nu

flatsamples = np.zeros((1, resolution))
flatsamplesparameters = []
i = 0
# flatsamples[i,:] = train_y.numpy()
flatsamples[i, :] = observed_pred_mean
flatsamplesparameters.append(test_x_m.numpy())
i += 1
maxindx = np.unravel_index(flatsamples.argmax(), flatsamples.shape)
max_parameters = flatsamplesparameters[maxindx[0]][maxindx[1]]
max_loglike = flatsamples.max()
maxpGBpredicted = scaletooriginal(max_parameters, boundaries_reduced)
if loglikelihood(maxpGBpredicted) > loglikelihood(maxpGB):
    maxpGB = maxpGBpredicted
    best_value = loglikelihood(maxpGB)
print(
    "pred",
    max_loglike,
    "true",
    loglikelihood(scaletooriginal(max_parameters, boundaries_reduced)),
    "max",
    loglikelihood(maxpGB),
    maxpGB,
)

start = time.time()
normalizer = sum(np.exp(observed_pred_mean-best_value))
flatsamples_normalized = np.exp(flatsamples[0]-best_value)/normalizer
mcmc_samples = np.zeros((resolution, len(parameters)))
mcmc_samples[0] = flatsamplesparameters[0][0]
previous_p = flatsamples_normalized[0]
for i in range(10**7-1):
    # print(flatsamples_normalized[i+1],previous_p,previous_p / flatsamples_normalized[i+1])
    if (flatsamples_normalized[i+1] / previous_p)**(1/2) > np.random.uniform():
        previous_p = flatsamples_normalized[i+1]
        mcmc_samples[i+1] = flatsamplesparameters[0][i+1]
    else:
        mcmc_samples[i+1] = mcmc_samples[i]
print('time MHMC', time.time()-start)
start = time.time()
# mcmc_samples_rescaled = np.zeros(np.shape(mcmc_samples))
# for k in range(len(mcmc_samples)):
#     rescaled = scaletooriginal(mcmc_samples[k], boundaries_reduced)
#     i = 0
#     for parameter in parameters:
#         mcmc_samples_rescaled[k,i]  =  rescaled[parameter]
#         i += 1
for parameter in parametersfd:
    i = 0
    if parameter in ["EclipticLatitude"]:
        mcmc_samples_rescaled[:,i] = np.arcsin((mcmc_samples[:,parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
    elif parameter in ["Inclination"]:
        mcmc_samples_rescaled[:,i] = np.arccos((mcmc_samples[:,parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
    else:
        mcmc_samples_rescaled[:,i] = (mcmc_samples[:,parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
        i += 1
print('time rescale', time.time()-start)
datS = np.zeros((resolution, 8))
datS[:,0] = mcmc_samples_rescaled[:,2]
datS[:,1] = mcmc_samples_rescaled[:,1]
datS[:,2] = mcmc_samples_rescaled[:,3]
datS[:,3] = np.log10(mcmc_samples_rescaled[:,4])
datS[:,4] = mcmc_samples_rescaled[:,5]
datS[:,5] = np.log10(mcmc_samples_rescaled[:,0])
datS[:,6] = mcmc_samples_rescaled[:,6]
datS[:,7] = mcmc_samples_rescaled[:,7]
fig2 =  corner.corner(datS,  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
                        color='#348ABD', use_math_test=True,\
                        levels=[0.9], title_kwargs={"fontsize": 12})
plt.show()

parameters_recorded[best_run]["Loglikelihood"].append(loglikelihood(maxpGB))
for parameter in parametersfd:
    parameters_recorded[best_run][parameter].append(maxpGB[parameter])
fig, ax = plt.subplots(2, 5, figsize=np.asarray(fig_size) * 2)
# plt.suptitle("Intermediate Parameters")
for n in range(len(parameters_recorded)):
    i = 0
    for parameter in parametersfd + ["Loglikelihood"]:
        j = 0
        if i > 3:
            j = 1
        if parameter in ["Frequency"]:
            ax[j, i % 5].plot(
                np.asarray(parameters_recorded[n][parameter]) * 10 ** 3,
                np.arange(0, len(parameters_recorded[n][parameter])),
            )
        else:
            ax[j, i % 5].plot(
                parameters_recorded[n][parameter],
                np.arange(0, len(parameters_recorded[n][parameter])),
            )
        i += 1
i = 0
for parameter in parametersfd + ["Log-likelihood"]:
    j = 0
    if i > 3:
        j = 1
    if parameter in ["Log-likelihood"]:
        ax[j, i % 5].axvline(x=loglikelihood(pGB), color="k", label="True")
    elif parameter in ["Frequency"]:
        ax[j, i % 5].axvline(x=pGB[parameter] * 10 ** 3, color="k", label="True")
    else:
        ax[j, i % 5].axvline(x=pGB[parameter], color="k", label="True")
    if parameter in ["Amplitude"]:
        ax[j, i % 5].set_xlabel(parameter, ha="right", va="top")
    elif parameter in ["Frequency"]:
        ax[j, i % 5].xaxis.set_major_locator(plt.MaxNLocator(3))
        ax[j, i % 5].set_xlabel(parameter + " (mHz)", ha="right", va="top")
    else:
        ax[j, i % 5].xaxis.set_major_locator(plt.MaxNLocator(3))
        ax[j, i % 5].set_xlabel(parameter)
    # if parameter in ['Log-likelihood']:
    #     ax[j,i%5].set_xlabel('$$')
    ax[j, i % 1].set_ylabel("Step")
    ax[j, i % 5].set_ylim(0, 100)
    if parameter in ["Log-likelihood"]:
        pass
    elif parameter == "EclipticLatitude":
        ax[j, i % 5].set_xlim(np.arcsin(boundaries[parameter][0]), np.arcsin(boundaries[parameter][1]))
    elif parameter == "Inclination":
        ax[j, i % 5].set_xlim(np.arccos(boundaries[parameter][1]), np.arccos(boundaries[parameter][0]))
    elif parameter in ["Frequency"]:
        ax[j, i % 5].set_xlim(boundaries[parameter][0] * 10 ** 3, boundaries[parameter][1] * 10 ** 3)
    else:
        ax[j, i % 5].set_xlim(boundaries[parameter][0], boundaries[parameter][1])
    i += 1
ax[0, 0].legend()
# plt.tight_layout()
plt.show()

Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator="synthlisa")
index_low = np.searchsorted(Xs.f, dataX.f[0])
Xs = Xs[index_low : index_low + len(dataX)]
Ys = Ys[index_low : index_low + len(dataY)]
Zs = Zs[index_low : index_low + len(dataZ)]

plt.figure(figsize=np.asarray(fig_size) * 2)
ax1 = plt.subplot(231)
# plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
ax1.plot(dataX.f * 1000, dataX.values.real, label="Data")
ax1.plot(Xs.f * 1000, Xs.values.real, label="GB")
ax2 = plt.subplot(232)
# plt.plot(dataY_training.f*1000,dataY_training.values, label='data')
ax2.plot(dataY.f * 1000, dataY.values.real, label="Data")
ax2.plot(Ys.f * 1000, Ys.values.real, label="True Response")
ax3 = plt.subplot(233)
# plt.plot(dataZ_training.f*1000,dataZ_training.values, label='data')
ax3.plot(dataZ.f * 1000, dataZ.values.real, label="Data")
ax3.plot(Zs.f * 1000, Zs.values.real, label="True Response")
ax4 = plt.subplot(234)
# plt.plot(dataX_training.f*1000,dataX_training.values.imag, label='data')
ax4.plot(dataX.f * 1000, dataX.values.imag, label="Data")
ax4.plot(Xs.f * 1000, Xs.values.imag, label="True Response")
ax5 = plt.subplot(235)
# plt.plot(dataY_training.f*1000,dataY_training.values.imag, label='data')
ax5.plot(dataY.f * 1000, dataY.values.imag, label="Data")
ax5.plot(Ys.f * 1000, Ys.values.imag, label="True Response")
ax6 = plt.subplot(236)
# plt.plot(dataZ_training.f*1000,dataZ_training.values.imag, label='data')
ax6.plot(dataZ.f * 1000, dataZ.values.imag, label="Data")
ax6.plot(Zs.f * 1000, Zs.values.imag, label="True Response")
Xs, Ys, Zs = GB.get_fd_tdixyz(template=maxpGB, oversample=4, simulator="synthlisa")
index_low = np.searchsorted(Xs.f, dataX.f[0])
Xs = Xs[index_low : index_low + len(dataX)]
Ys = Ys[index_low : index_low + len(dataY)]
Zs = Zs[index_low : index_low + len(dataZ)]
ax1.plot(Xs.f * 1000, Xs, label="Found GB")
ax1.xaxis.set_major_locator(plt.MaxNLocator(3))
ax1.set_xlabel("f (mHz)")
ax1.set_ylabel("X-TDI real (1/Hz)")
ax1.legend(loc="upper left", numpoints=0.5)
ax2.plot(Ys.f * 1000, Ys, label="Found Response")
ax2.xaxis.set_major_locator(plt.MaxNLocator(3))
ax2.set_xlabel("f (mHz)")
ax2.set_ylabel("Y-TDI real (1/Hz)")
# ax2.legend()
ax3.plot(Zs.f * 1000, Zs, label="Found Response")
ax3.xaxis.set_major_locator(plt.MaxNLocator(3))
ax3.set_xlabel("f (mHz)")
ax3.set_ylabel("Z-TDI real (1/Hz)")
# ax3.legend()
ax4.plot(Xs.f * 1000, Xs.imag, label="Found Response")
ax4.xaxis.set_major_locator(plt.MaxNLocator(3))
ax4.set_xlabel("f (mHz)")
ax4.set_ylabel("X-TDI imag (1/Hz)")
# ax4.legend()
ax5.plot(Ys.f * 1000, Ys.imag, label="Found Response")
ax5.xaxis.set_major_locator(plt.MaxNLocator(3))
ax5.set_xlabel("f (mHz)")
ax5.set_ylabel("Y-TDI imag (1/Hz)")
# ax5.legend()
ax6.plot(Zs.f * 1000, Zs.imag, label="Found Response")
ax6.xaxis.set_major_locator(plt.MaxNLocator(3))
ax6.set_xlabel("f (mHz)")
ax6.set_ylabel("Z-TDI imag (1/Hz)")
# ax6.legend()
plt.tight_layout()
plt.show()

# Get into evaluation (predictive posterior) mode
model.eval()
likelihood.eval()
observed_pred = {}

number_of_test_samples = 100
test_x = {}
test_y = {}
for parameter in parametersfd:
    test_samples = sampler(
        number_of_test_samples,
        parameters,
        maxpGB,
        boundaries_reduced,
        p1,
        uniform=True,
        only=True,
        onlyparameter=parameter,
    )
    test_x[parameter] = np.zeros((number_of_test_samples, len(parameters)))
    i = 0
    for name in parametersfd:
        test_x[parameter][:, i] = test_samples[name]
        i += 1
    test_y[parameter] = test_samples["Likelihood"]
    test_x[parameter] = torch.from_numpy(test_x[parameter]).float()
    test_x[parameter] = test_x[parameter]
    test_y[parameter] = torch.from_numpy(test_y[parameter]).float()
number_of_test_samples = 500
test_x2 = {}
test_y2 = {}

for parameter in parametersfd:
    test_samples = sampler(
        number_of_test_samples,
        parameters,
        maxpGB,
        boundaries,
        p1,
        uniform=True,
        only=True,
        onlyparameter=parameter,
    )
    test_x2[parameter] = np.zeros((number_of_test_samples, len(parameters)))
    i = 0
    for name in parametersfd:
        test_x2[parameter][:, i] = test_samples[name]
        i += 1
    test_y2[parameter] = test_samples["Likelihood"]
    test_x2[parameter] = torch.from_numpy(test_x2[parameter]).float()
    test_x2[parameter] = test_x2[parameter]
    test_y2[parameter] = torch.from_numpy(test_y2[parameter]).float()
number_of_test_samples2d = 50 ** 2

parameter = "EclipticLongitude"
parameters_reduced = [
    "EclipticLatitude",
    "Frequency",
    "Inclination",
    "Amplitude",
    'InitialPhase',
    'Polarization'
]
observed_pred_mean = {}
pred_sigma = {}
for n in range(len(parameters_reduced)):
    for parameter2 in parameters_reduced:
        print(parameter+parameter2)
        if parameter != parameter2:
            test_samples = sampler(
                number_of_test_samples2d,
                parameters,
                maxpGB,
                boundaries_reduced,
                p1,
                uniform=True,
                twoD=True,
                onlyparameter=parameter,
                secondparameter=parameter2,
                calculate_loglikelihood=False,
            )
            test_x[parameter + parameter2] = np.zeros((number_of_test_samples2d, len(parameters)))
            i = 0
            for name in parametersfd:
                test_x[parameter + parameter2][:, i] = test_samples[name]
                i += 1
            test_x[parameter + parameter2] = torch.from_numpy(test_x[parameter + parameter2]).float()

            observed_pred_sk, pred_sigma[parameter + parameter2] = gpr.predict(test_x[parameter + parameter2],return_std=True)
            observed_pred_mean[parameter + parameter2] = observed_pred_sk*sigma +nu
            pred_sigma[parameter + parameter2] = pred_sigma[parameter + parameter2]*sigma +nu
            # with torch.no_grad(), gpytorch.settings.fast_pred_var():
            #     observed_pred[parameter + parameter2] = likelihood(model(test_x[parameter + parameter2]))
    parameter = parameters_reduced[0]
    parameters_reduced.pop(0)


for parameter in parametersfd:
    # Make predictions by feeding model through likelihood
    observed_pred_sk, pred_sigma[parameter] = gpr.predict(test_x[parameter],return_std=True)
    observed_pred_mean[parameter] = observed_pred_sk*sigma +nu
    pred_sigma[parameter] = pred_sigma[parameter]*sigma +nu
    # with torch.no_grad(), gpytorch.settings.fast_pred_var():
    #     observed_pred[parameter] = likelihood(model(test_x[parameter]))
    #     observed_pred_print = (observed_pred[parameter].mean.numpy() * sigma) + nu

    print(
        "sqrt(MSE) ",
        parameter,
        np.sqrt(mean_squared_error(test_y[parameter].numpy(), observed_pred_mean[parameter])),
    )
    # if parameter != parameter2:
    #     with torch.no_grad(), gpytorch.settings.fast_pred_var():
    #         observed_pred[parameter + parameter2] = likelihood(model(test_x[parameter + parameter2]))

# parameter = 'random'
# with torch.no_grad(), gpytorch.settings.fast_pred_var():
#     observed_pred[parameter] = likelihood(model(test_x[parameter]))
# observed_pred_print = (observed_pred[parameter]*sigma)+nu
# print('sqrt(MSE) ','random', np.sqrt(mean_squared_error(test_y[parameter].numpy(),observed_pred_print.mean.cpu().numpy())))


for parameter in parameters:
    if parameter in ["EclipticLatitude"]:
        samples[parameter] = np.arcsin((samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        test_samples[parameter] = np.arcsin((test_samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
    elif parameter in ["Inclination"]:
        samples[parameter] = np.arccos((samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        test_samples[parameter] = np.arccos((test_samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
    else:
        samples[parameter] = (samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
        test_samples[parameter] = (test_samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]


train_x = train_x.cpu()
train_y = train_y.cpu()
# train_y = (train_y*sigma)+nu

pGB01 = {}
for parameter in parameters:
    if parameter in ["EclipticLatitude"]:
        pGB01[parameter] = (np.sin(pGB[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
    elif parameter in ["Inclination"]:
        pGB01[parameter] = (np.cos(pGB[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
    else:
        pGB01[parameter] = (pGB[parameter] - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])


fig, ax = plt.subplots(2, 4, figsize=np.asarray(fig_size) * 2.5)
# plt.suptitle("loglikelihood")
i = 0
mean = {}
for parameter in parametersfd:
    j = 0
    if i > 3:
        j = 1
    test_x_rescaled = []
    for k in range(len(test_x[parameter])):
        test_x_rescaled.append(scaletooriginal(test_x[parameter][k].numpy(), boundaries_reduced)[parameter])
    test_x2_rescaled = []
    for k in range(len(test_x2[parameter])):
        test_x2_rescaled.append(scaletooriginal(test_x2[parameter][k].numpy(), boundaries)[parameter])
    # Plot training data as black stars
    if parameter == "Frequency":
        ax[j, i % 4].axvline(x=pGB[parameter] * 10 ** 3, color="k", label="GB")
        ax[j, i % 4].plot(
            np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            np.exp(test_y[parameter].numpy()[1:]-best_value)/normalizer,
            "g",linewidth=3, alpha=0.5,
            label="True",
        )
        ax[j, i % 4].plot(
            np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            np.exp(observed_pred_mean[parameter][1:]-best_value)/normalizer,
            "b",
            label="Mean",
        )
        # ax[j, i % 4].fill_between(
        #     np.asarray(test_x_rescaled[1:]) * 10 ** 3,
        #     np.exp(observed_pred_mean[parameter][1:]+pred_sigma[parameter][1:]-best_value)/normalizer,
        #     np.exp(observed_pred_mean[parameter][1:]-pred_sigma[parameter][1:]-best_value)/normalizer,
        #     alpha=0.5,
        #     label="Confidence",
        # )
    else:
        ax[j, i % 4].axvline(x=pGB[parameter], color="k", label="GB")
        ax[j, i % 4].plot(test_x_rescaled[1:], np.exp(test_y[parameter].numpy()[1:]-best_value)/normalizer, "g",linewidth=5, alpha=0.5, label="True")
        ax[j, i % 4].plot(test_x_rescaled[1:], np.exp(observed_pred_mean[parameter][1:]-best_value)/normalizer, "b", label="Mean")
        # ax[j, i % 4].fill_between(
        #     test_x_rescaled[1:],
        #     np.exp(observed_pred_mean[parameter][1:]+pred_sigma[parameter][1:]-best_value)/normalizer,
        #     np.exp(observed_pred_mean[parameter][1:]-pred_sigma[parameter][1:]-best_value)/normalizer,
        #     alpha=0.5,
        #     label="Confidence",
        # )
    ax[0, 0].legend()
    if parameter in ["Amplitude"]:
        ax[j, i % 4].set_xlabel(parameter, ha="right", va="top")
    elif parameter in ["Frequency"]:
        ax[j, i % 4].set_xlabel(parameter + " (mHz)", ha="right", va="top")
    else:
        ax[j, i % 4].set_xlabel(parameter)
    ax[j, i % 4].set_ylabel("Likelihood")
    if parameter == "EclipticLatitude":
        ax[j, i % 4].set_xlim(np.arcsin(boundaries_reduced[parameter][0]), np.arcsin(boundaries_reduced[parameter][1]))
    elif parameter == "Inclination":
        ax[j, i % 4].set_xlim(np.arccos(boundaries_reduced[parameter][1]), np.arccos(boundaries_reduced[parameter][0]))
    elif parameter in ["Frequency"]:
        ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0] * 10 ** 3, boundaries_reduced[parameter][1] * 10 ** 3)
    else:
        ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
    i += 1
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(2, 4, figsize=np.asarray(fig_size) * 2.5)
# plt.suptitle("loglikelihood")
i = 0
for parameter in parametersfd:
    j = 0
    if i > 3:
        j = 1
    with torch.no_grad():
        # Get upper and lower confidence bounds
        test_x_rescaled = []
        for k in range(len(test_x[parameter])):
            test_x_rescaled.append(scaletooriginal(test_x[parameter][k].numpy(), boundaries_reduced)[parameter])
        test_x2_rescaled = []
        for k in range(len(test_x2[parameter])):
            test_x2_rescaled.append(scaletooriginal(test_x2[parameter][k].numpy(), boundaries)[parameter])
        # Plot training data as black stars
        if parameter == "Frequency":
            ax[j, i % 4].axvline(x=pGB[parameter] * 10 ** 3, color="k", label="GB")
            ax[j, i % 4].plot(
                np.asarray(test_x2_rescaled[1:]) * 10 ** 3,
                test_y2[parameter].numpy()[1:],
                "g",
                label="True",
            )
            ax[j, i % 4].plot(
                np.asarray(test_x_rescaled[1:]) * 10 ** 3,
                observed_pred_mean[parameter][1:],
                "b",
                label="Mean",
            )
            # ax[j, i % 4].fill_between(
            #     np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            #     observed_pred_mean[parameter][1:]+pred_sigma[parameter][1:],
            #     observed_pred_mean[parameter][1:]-pred_sigma[parameter][1:],
            #     alpha=0.5,
            #     label="Confidence",
            # )
        else:
            ax[j, i % 4].axvline(x=pGB[parameter], color="k", label="GB")
            ax[j, i % 4].plot(test_x2_rescaled[1:], test_y2[parameter].numpy()[1:], "g", label="True")
            ax[j, i % 4].plot(test_x_rescaled[1:], observed_pred_mean[parameter][1:], "b", label="Mean")
            # ax[j, i % 4].fill_between(
            #     test_x_rescaled[1:],
            #     observed_pred_mean[parameter][1:]+pred_sigma[parameter][1:],
            #     observed_pred_mean[parameter][1:]-pred_sigma[parameter][1:],
            #     alpha=0.5,
            #     label="Confidence",
            # )
        ax[0, 0].legend()
    if parameter in ["Amplitude"]:
        ax[j, i % 4].set_xlabel(parameter, ha="right", va="top")
    elif parameter in ["Frequency"]:
        ax[j, i % 4].set_xlabel(parameter + " (mHz)", ha="right", va="top")
    else:
        ax[j, i % 4].set_xlabel(parameter)
    ax[j, i % 4].set_ylabel("Log-likelihood")
    # if parameter == "EclipticLatitude":
    #     ax[j, i % 4].set_xlim(np.arcsin(boundaries_reduced[parameter][0]), np.arcsin(boundaries_reduced[parameter][1]))
    # elif parameter == "Inclination":
    #     ax[j, i % 4].set_xlim(np.arccos(boundaries_reduced[parameter][1]), np.arccos(boundaries_reduced[parameter][0]))
    # elif parameter in ["Frequency"]:
    #     ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0] * 10 ** 3, boundaries_reduced[parameter][1] * 10 ** 3)
    # else:
    #     ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
    i += 1
plt.tight_layout()
plt.show()
# fig, ax = plt.subplots(2, 4,figsize=(15,15))
# plt.suptitle("loglikelihood true")
# i = 0
# for parameter in parametersfd:
#     if parameter != parameter2:
#         j = 0
#         if i > 3:
#             j = 1
#         with torch.no_grad():
#             ax[j,i%4].axvline(x=pGB01[parameter], color='k')
#             ax[j,i%4].axhline(y=pGB01[parameter2], color='k')
#             im = ax[j,i%4].scatter(test_x[parameter+parameter2].numpy()[:,i],test_x[parameter+parameter2].numpy()[:,parametersfd.index(parameter2)],c=test_y[parameter+parameter2].numpy()[:])
#             ax[j,i%4].set_xlabel(parameter)
#             ax[j,i%4].set_ylabel(parameter2)
#             fig.colorbar(im, ax=ax[j,i%4])
#     else:
#         ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], test_y[parameter].numpy()[1:], 'g.')
#         ax[j,i%4].set_xlabel(parameter)
#         ax[j,i%4].set_ylabel('loglikelihood')
#         ax[j,i%4].legend(['True'])
#     i += 1

fig, ax = plt.subplots(2, 4, figsize=(15, 15))
plt.suptitle("loglikelihood predicted mean")
i = 0
for parameter in parametersfd:
    if parameter != parameter2:
        j = 0
        if i > 3:
            j = 1
        with torch.no_grad():
            # Get upper and lower confidence bounds
            mean2 = observed_pred[parameter + parameter2].mean
            mean[parameter + parameter2] = (mean2 * sigma) + nu
            test_x_rescaled1 = []
            test_x_rescaled2 = []
            for k in range(len(test_x[parameter + parameter2])):
                test_x_rescaled1.append(scaletooriginal(test_x[parameter + parameter2][k].numpy(), boundaries_reduced)[parameter])
                test_x_rescaled2.append(scaletooriginal(test_x[parameter + parameter2][k].numpy(), boundaries_reduced)[parameter2])
            ax[j, i % 4].axvline(x=pGB[parameter], color="k")
            ax[j, i % 4].axhline(y=pGB[parameter2], color="k")
            im = ax[j, i % 4].scatter(test_x_rescaled1, test_x_rescaled2, c=mean[parameter + parameter2][:])
            ax[j, i % 4].set_xlabel(parameter)
            ax[j, i % 4].set_ylabel(parameter2)
            fig.colorbar(im, ax=ax[j, i % 4])

            if parameter == "EclipticLatitude":
                ax[j, i % 4].set_xlim(
                    np.arcsin(boundaries_reduced[parameter][0]),
                    np.arcsin(boundaries_reduced[parameter][1]),
                )
            elif parameter == "Inclination":
                ax[j, i % 4].set_xlim(
                    np.arccos(boundaries_reduced[parameter][1]),
                    np.arccos(boundaries_reduced[parameter][0]),
                )
            elif parameter in ["Frequency"]:
                ax[j, i % 4].set_xlim(
                    boundaries_reduced[parameter][0] * 10 ** 3,
                    boundaries_reduced[parameter][1] * 10 ** 3,
                )
            else:
                ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
            if parameter2 == "EclipticLatitude":
                ax[j, i % 4].set_ylim(
                    np.arcsin(boundaries_reduced[parameter2][0]),
                    np.arcsin(boundaries_reduced[parameter2][1]),
                )
            elif parameter2 == "Inclination":
                ax[j, i % 4].set_ylim(
                    np.arccos(boundaries_reduced[parameter2][1]),
                    np.arccos(boundaries_reduced[parameter2][0]),
                )
            elif parameter2 in ["Frequency"]:
                ax[j, i % 4].set_ylim(
                    boundaries_reduced[parameter2][0] * 10 ** 3,
                    boundaries_reduced[parameter2][1] * 10 ** 3,
                )
            else:
                ax[j, i % 4].set_ylim(boundaries_reduced[parameter2][0], boundaries_reduced[parameter2][1])
    else:
        with torch.no_grad():
            lower, upper = observed_pred[parameter].confidence_region()
            lower = lower.cpu()
            upper = upper.cpu()
            lower = (lower * sigma) + nu
            upper = (upper * sigma) + nu
            ax[j, i % 4].plot(test_x[parameter].numpy()[1:, i], test_y[parameter].numpy()[1:], "g.")
            ax[j, i % 4].plot(test_x[parameter].numpy()[1:, i], mean[parameter].numpy()[1:], "b.")
            ax[j, i % 4].fill_between(
                test_x[parameter].numpy()[1:, i],
                lower.numpy()[1:],
                upper.numpy()[1:],
                alpha=0.5,
            )
            ax[j, i % 4].set_xlabel(parameter)
            ax[j, i % 4].set_ylabel("loglikelihood")
            ax[j, i % 4].legend(["True", "Mean", "Confidence"])
    i += 1
plt.show()


fig, ax = plt.subplots(7, 7, figsize=np.asarray(fig_size)*3)
ax[0,1].axis('off')
ax[0,2].axis('off')
ax[0,3].axis('off')
ax[0,4].axis('off')
ax[1,2].axis('off')
ax[1,3].axis('off')
ax[1,4].axis('off')
ax[2,3].axis('off')
ax[2,4].axis('off')
ax[3,4].axis('off')
i = 0
parameter = "EclipticLatitude"
parameters_reduced = [
    "EclipticLongitude",
    "EclipticLatitude",
    "Frequency",
    "Inclination",
    "Amplitude",
    'InitialPhase',
    'Polarization'
]
for parameter in parameters_reduced:
    test_x_rescaled = []
    for k in range(len(test_x[parameter])):
        test_x_rescaled.append(scaletooriginal(test_x[parameter][k].numpy(), boundaries_reduced)[parameter])
    test_x2_rescaled = []
    for k in range(len(test_x2[parameter])):
        test_x2_rescaled.append(scaletooriginal(test_x2[parameter][k].numpy(), boundaries)[parameter])
    # Plot training data as black stars
    if parameter == "Frequency":
        ax[i,i].axvline(x=pGB[parameter] * 10 ** 3, color="k", label="GB")
        ax[i,i].plot(
            np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            np.exp(test_y[parameter].numpy()[1:]-best_value)/normalizer,
            "g",
            label="Target",
        )
        ax[i,i].plot(
            np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            np.exp(observed_pred_mean[parameter][1:]-best_value)/normalizer,
            "b",
            label="Mean",
        )
        # ax[i,i].fill_between(
        #     np.asarray(test_x_rescaled[1:]) * 10 ** 3,
        #     np.exp(lower.numpy()[1:]-best_value)/normalizer,
        #     np.exp(upper.numpy()[1:]-best_value)/normalizer,
        #     alpha=0.5,
        #     label="Confidence",
        # )
    else:
        ax[i,i].axvline(x=pGB[parameter], color="r", label="GB")
        ax[i,i].plot(test_x_rescaled[1:], np.exp(test_y[parameter].numpy()[1:]-best_value)/normalizer, "g", label="Target")
        ax[i,i].plot(test_x_rescaled[1:], np.exp(observed_pred_mean[parameter][1:]-best_value)/normalizer, "b", label="Mean")
        # ax[i,i].fill_between(
        #     test_x_rescaled[1:],
        #     np.exp(lower.numpy()[1:]-best_value)/normalizer,
        #     np.exp(upper.numpy()[1:]-best_value)/normalizer,
        #     alpha=0.5,
        #     label="Confidence",
        # )
    ax[4, 4].legend()
    ax[4,4].set_xlabel("Amplitude", ha="right", va="top")
    ax[0,0].set_ylabel('Likelihood')
    if i != 4:
        ax[i,i].xaxis.set_ticklabels([])
    if i != 0:
        ax[i,i].yaxis.set_ticklabels([])
    # ax[i,i].set_ylabel("Likelihood")
    if parameter == "EclipticLatitude":
        ax[i,i].set_xlim(np.arcsin(boundaries_reduced[parameter][0]), np.arcsin(boundaries_reduced[parameter][1]))
    elif parameter == "Inclination":
        ax[i,i].set_xlim(np.arccos(boundaries_reduced[parameter][1]), np.arccos(boundaries_reduced[parameter][0]))
    elif parameter in ["Frequency"]:
        ax[i,i].set_xlim(boundaries_reduced[parameter][0] * 10 ** 3, boundaries_reduced[parameter][1] * 10 ** 3)
    else:
        ax[i,i].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
    i += 1
parameters_reduced = [
    "EclipticLatitude",
    "Frequency",
    "Inclination",
    "Amplitude",
    'InitialPhase',
    'Polarization'
]
parameter = "EclipticLongitude"
for n in range(len(parameters_reduced)):
    print(n)
    print(parameter2)
    i = n+1
    for parameter2 in parameters_reduced:
        test_x_rescaled1 = []
        test_x_rescaled2 = []
        for k in range(len(test_x[parameter + parameter2])):
            test_x_rescaled1.append(scaletooriginal(test_x[parameter + parameter2][k].numpy(), boundaries_reduced)[parameter])
            test_x_rescaled2.append(scaletooriginal(test_x[parameter + parameter2][k].numpy(), boundaries_reduced)[parameter2])
        ax[i,n].axvline(x=pGB[parameter], color="r")
        ax[i,n].axhline(y=pGB[parameter2], color="r")
        ax[i,n].scatter(test_x_rescaled1, test_x_rescaled2, c=np.exp(observed_pred_mean[parameter + parameter2][:]-best_value)/normalizer)
        if parameter == 'Frequency':
            ax[i,n].scatter(np.asarray(test_x_rescaled1)*1e3, test_x_rescaled2, c=np.exp(observed_pred_mean[parameter + parameter2][:]-best_value)/normalizer)
            ax[i,n].axvline(x=pGB[parameter]*1e3, color="r")
        if parameter2 == 'Frequency':
            ax[i,n].scatter(test_x_rescaled1, np.asarray(test_x_rescaled2)*1e3, c=np.exp(observed_pred_mean[parameter + parameter2][:]-best_value)/normalizer)
        ax[i,n].axhline(y=pGB[parameter2]*1e3, color="r")
        if i == 6:
            ax[i,n].set_xlabel(parameter)
            if parameter in ["Frequency"]:
                ax[i,n].set_xlabel(parameter + " (mHz)", ha="right", va="top" )
                ax[i,n].xaxis.set_major_locator(plt.MaxNLocator(2))
        else:
            ax[i,n].xaxis.set_ticklabels([])
        if n == 0:
            ax[i,n].set_ylabel(parameter2)
            if parameter in ["Frequency"]:
                ax[i,n].set_ylabel(parameter + " (mHz)")
        else:
            ax[i,n].yaxis.set_ticklabels([])
        # fig.colorbar(im, ax=ax[i,n])

        if parameter == "EclipticLatitude":
            ax[i,n].set_xlim(
                np.arcsin(boundaries_reduced[parameter][0]),
                np.arcsin(boundaries_reduced[parameter][1]),
            )
        elif parameter == "Inclination":
            ax[i,n].set_xlim(
                np.arccos(boundaries_reduced[parameter][1]),
                np.arccos(boundaries_reduced[parameter][0]),
            )
        elif parameter in ["Frequency"]:
            ax[i,n].set_xlim(
                boundaries_reduced[parameter][0] * 10 ** 3,
                boundaries_reduced[parameter][1] * 10 ** 3,
            )
        else:
            ax[i,n].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
        if parameter2 == "EclipticLatitude":
            ax[i,n].set_ylim(
                np.arcsin(boundaries_reduced[parameter2][0]),
                np.arcsin(boundaries_reduced[parameter2][1]),
            )
        elif parameter2 == "Inclination":
            ax[i,n].set_ylim(
                np.arccos(boundaries_reduced[parameter2][1]),
                np.arccos(boundaries_reduced[parameter2][0]),
            )
        elif parameter2 in ["Frequency"]:
            ax[i,n].set_ylim(
                boundaries_reduced[parameter2][0] * 10 ** 3,
                boundaries_reduced[parameter2][1] * 10 ** 3,
            )
        else:
            ax[i,n].set_ylim(boundaries_reduced[parameter2][0], boundaries_reduced[parameter2][1])
        i += 1
    parameter = parameters_reduced[0]
    parameters_reduced.pop(0)
plt.show()

prediction = observed_pred["Frequency"].mean.cpu().numpy()
n_bin = 50
i = 0
fig, axes = plt.subplots(2, 4, figsize=(15, 15))
plt.suptitle("sampled posterior")
for parameter in parameters:
    j = 0
    if i > 3:
        j = 1
    axes[j, i % 4].axvline(x=pGB[parameter], color="r")
    axes[j, i % 4].plot(samples[parameter], samples["Likelihood"], ".", label="train")
    axes[j, i % 4].plot(test_samples[parameter], test_samples["Likelihood"], ".", label="true")
    axes[j, i % 4].errorbar(
        test_samples[parameter][1:],
        prediction[1:],
        prediction[1:] - lower.numpy()[1:],
        fmt=".",
        label="prediction",
    )
    # plt.fill_between(test_samples[parameter][1:],prediction[1:]-p_std[1:],prediction[1:]+p_std[1:], color='g', alpha= 0.4)
    if i == 0:
        axes[j, i % 4].set_ylabel("log-Likelihood")
    # plt.subplot(3,number_of_parameters,i+number_of_parameters)
    # plt.axvline(x=pGB[parameter], color='r')
    # n, bins, patches = plt.hist(samples[parameter], n_bin, density=True, facecolor='k', alpha=0.5)

    # plt.subplot(3,number_of_parameters,i+2*number_of_parameters)
    # if parameter == 'Frequency':
    #     plt.axvline(x=pGB[parameter]*1000, color='r')
    #     plt.plot(samples[parameter]*1000, range(number_of_samples), 'k')
    # else:
    #     plt.axvline(x=pGB[parameter], color='r')
    #     plt.plot(samples[parameter], range(number_of_samples), 'k')
    axes[j, i % 4].set_xlabel(parameter)
    i += 1
plt.legend()

name = "Frequency"
plt.figure(figsize=(10, 8))
plt.suptitle("sampled posterior")
plt.subplot(1, 1, 1)
plt.axvline(x=pGB[name], color="r")
plt.plot(samples[name], samples["Likelihood"], ".", label="train", zorder=1)
plt.plot(test_samples[name], test_samples["Likelihood"], ".", label="true")
plt.plot(test_samples[name][1:], prediction[1:], label="prediction")
plt.fill_between(test_samples[name][1:], lower.numpy()[1:], upper.numpy()[1:], color="g", alpha=0.4)
plt.ylabel("log-Likelihood")
plt.xlabel(name)
plt.legend()
# plt.figure()
# plt.scatter(samples['Amplitude'],samples['Frequency'],c=samples['Likelihood'])
# plt.xlim(max(samples['Amplitude']),max(samples['Amplitude']))
# plt.ylim(max(samples['Frequency']),max(samples['Frequency']))
# plt.colorbar()
# plt.figure()
# # plt.title('sqrt(MSE)',np.round(np.sqrt(mean_squared_error(test_y,prediction)),2))
# plt.scatter(samples['EclipticLatitude'],samples['EclipticLongitude'],c=samples['Likelihood'])
# plt.xlim(max(samples['EclipticLatitude']),max(samples['EclipticLatitude']))
# plt.ylim(max(samples['EclipticLongitude']),max(samples['EclipticLongitude']))
# plt.colorbar()
plt.show()
