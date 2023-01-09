#%%
from re import A
# from matplotlib.lines import _LineStyle
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy
from scipy.optimize import differential_evolution
import numpy as np
import xarray as xr
import time
from copy import deepcopy
import multiprocessing as mp
import pandas as pd
import os
import h5py
import sys
sys.path.append('/cluster/home/sstrub/Repositories/LDC/lib/lib64/python3.8/site-packages/ldc-0.1-py3.8-linux-x86_64.egg')

from astropy import units as u
import astropy.coordinates as coord

from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries, window
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr





# get current directory
path = os.getcwd()
 
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)
Radler = True
if Radler:
    DATAPATH = grandparent+"/LDC/Radler/data/"
    SAVEPATH = grandparent+"/LDC/pictures/LDC1-4/"
    SAVEPATH = grandparent+"/LDC/Radler/LDC1-4_evaluation/"
else:
    DATAPATH = grandparent+"/LDC/Sangria/data/"
    SAVEPATH = grandparent+"/LDC/pictures/Sangria/"


# save_name = 'Sangria_1_full_cut'
save_name = 'LDC1-4_2_optimized_second'
# save_name = 'LDC1-4_half_year'
# save_name2 = 'LDC1-4_half_year'

# duration = '3932160'
# duration = '7864320'
# duration = '15728640'
duration = '31457280'
save_name2 = 'Montana2022_'+duration

Tobs = int(duration)

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

end_string = 'amplitude_considered_sqrt'
# end_string = 'correlation'
found_sources_matched_flat_df = pd.read_pickle(SAVEPATH+'/found_sources_matched' +save_name+end_string+'_df')
found_sources_not_matched_flat_df = pd.read_pickle(SAVEPATH+'/found_sources_not_matched' +save_name+end_string+'_df')
pGB_injected_matched_flat_df = pd.read_pickle(SAVEPATH+'/injected_matched_windows' +save_name+end_string+'_df')
pGB_injected_not_matched_flat_df = pd.read_pickle(SAVEPATH+'/injected_not_matched_windows' +save_name+end_string+'_df')
correlation_list = np.load(SAVEPATH+'correlation_list' +save_name+end_string+'.npy', allow_pickle=True)

number_of_found_signals_not_matched = len(found_sources_not_matched_flat_df)
number_of_matched_signals = len(found_sources_matched_flat_df)
number_of_found_signals = number_of_matched_signals + number_of_found_signals_not_matched
number_of_injected_signals = len(pGB_injected_matched_flat_df) + len(pGB_injected_not_matched_flat_df)

print(number_of_matched_signals ,'matched signals out of', number_of_found_signals, 'found signals')
print('matched signals/found signals:', np.round(number_of_matched_signals/number_of_found_signals,2))

found_sources_matched_flat_df2 = pd.read_pickle(SAVEPATH+'/found_sources_matched' +save_name2+end_string+'_df')
found_sources_not_matched_flat_df2 = pd.read_pickle(SAVEPATH+'/found_sources_not_matched' +save_name2+end_string+'_df')
pGB_injected_matched_flat_df2 = pd.read_pickle(SAVEPATH+'/injected_matched_windows' +save_name2+end_string+'_df')
correlation_list2 = np.load(SAVEPATH+'correlation_list' +save_name2+end_string+'.npy', allow_pickle=True)

number_of_found_signals_not_matched2 = len(found_sources_not_matched_flat_df2)
number_of_matched_signals2 = len(found_sources_matched_flat_df2)
number_of_found_signals2 = number_of_matched_signals2 + number_of_found_signals_not_matched2

print(number_of_matched_signals2 ,'matched signals out of', number_of_found_signals2, 'found signals 2')
print('matched signals/found signals:', np.round(number_of_matched_signals2/number_of_found_signals2,2))

#### plot SNR - frequency
markersize = 3
alpha = 0.5
parameter_to_plot = 'IntrinsicSNR'
fig = plt.figure()
# plt.plot(pGB_injected_flat_df['Frequency']*10**3,pGB_injected_flat_df['IntrinsicSNR'], '.', color= colors[0], label = 'Injected', markersize= markersize, alpha = alpha)
# plt.plot(pGB_injected_matched_flat_df['Frequency']*10**3,pGB_injected_matched_flat_df['IntrinsicSNR'], '.', color= colors[1], label = 'Injected matched', markersize= markersize, alpha = alpha)
# plt.plot(pGB_injected_flat_df_high_SNR['Frequency']*10**3,pGB_injected_flat_df_high_SNR['IntrinsicSNR'],'.', color= colors[1], markersize= markersize, label = 'Injected SNR > 10', alpha = alpha)
plt.plot(found_sources_matched_flat_df['Frequency']*10**3,found_sources_matched_flat_df['IntrinsicSNR'],'g.', label = 'Found matched', markersize= markersize, alpha = alpha)
# plt.plot(pGB_injected_not_matched_flat_df['Frequency']*10**3,pGB_injected_not_matched_flat_df['IntrinsicSNR'], '+', color = 'r', label = 'Injected not matched', markersize= markersize, alpha = alpha)
plt.plot(found_sources_not_matched_flat_df['Frequency']*10**3,found_sources_not_matched_flat_df['IntrinsicSNR'],'o',color= 'blue',  markerfacecolor='None', markersize= markersize, label = 'Found not matched', alpha = alpha)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('f (mHz)')
plt.xlim(0.3,100)
plt.ylim(0.1,2000)
if parameter_to_plot == 'IntrinsicSNR':
    plt.ylabel('Intrinsic SNR')
else:
    plt.ylabel(parameter_to_plot)    
plt.legend(markerscale=4, loc = 'lower right')
plt.savefig(SAVEPATH+'/Evaluation/'+parameter_to_plot+save_name+'injected_not_matched_found_matched_found_not_matched'+end_string,dpi=300,bbox_inches='tight')
plt.show()

def get_errors(pGB_injected_matched_flat_df, found_sources_matched_flat_df):
    error = {}
    for parameter in parameters:
        error[parameter] = []
        # found_sources_matched_flat_df[parameter+'Error'] = []
        for i in range(len(pGB_injected_matched_flat_df)):
            if parameter == 'EclipticLongitude':
                if pGB_injected_matched_flat_df[parameter][i] > np.pi:
                    pGB_injected_matched_flat_df[parameter][i] -= 2*np.pi
            if parameter in ['EclipticLongitude', 'EclipticLatitude', 'Inclination', 'InitialPhase', 'Polarization']:
                error[parameter].append(np.abs(np.arcsin(np.sin(pGB_injected_matched_flat_df[parameter][i] - found_sources_matched_flat_df[parameter][i]))))
            elif parameter == 'Amplitude':
                # error[parameter].append(np.log10(pGB_injected_matched_flat_df[parameter][i]) - np.log10(found_sources_matched_flat_df[parameter][i]))
                error[parameter].append(np.abs(pGB_injected_matched_flat_df[parameter][i] - found_sources_matched_flat_df[parameter][i])/pGB_injected_matched_flat_df[parameter][i])
            # found_sources_matched_flat_df[parameter+'Error'].append(error[parameter][i])
            else:
                error[parameter].append(pGB_injected_matched_flat_df[parameter][i] - found_sources_matched_flat_df[parameter][i])
        found_sources_matched_flat_df[parameter+'Error'] = error[parameter]
    return found_sources_matched_flat_df

found_sources_matched_flat_df = get_errors(pGB_injected_matched_flat_df, found_sources_matched_flat_df)

found_sources_matched_flat_df2 = get_errors(pGB_injected_matched_flat_df2, found_sources_matched_flat_df2)

##### plot SNR to error
# for parameter in  parameters +['Skylocation']:
# for parameter in  ['SkylocationError']:
x_parameter = 'Frequency'
# x_parameter = 'IntrinsicSNR'
for parameter in  ['Frequency']:
    # x_parameter = parameter
    markersize = 4
    alpha = 0.5
    fig = plt.figure()
    # plt.plot(found_sources_matched_flat_df[parameter+'Error']/found_sources_matched_flat_df[parameter],found_sources_matched_flat_df[y_parameter],'.', label = 'found', markersize=markersize)
    # if parameter in ['Frequency', 'FrequencyDerivative']:
    #     plt.scatter(found_sources_matched_flat_df[x_parameter],found_sources_matched_flat_df[parameter+'Error']/found_sources_matched_flat_df['Frequency'], c=found_sources_matched_flat_df['Frequency'])
    # else:
    plt.plot(found_sources_matched_flat_df[x_parameter],found_sources_matched_flat_df[parameter+'Error'],'.', alpha =0.4, markersize=6)
        # plt.scatter(found_sources_matched_flat_df[x_parameter],found_sources_matched_flat_df[parameter], c=found_sources_matched_flat_df['Frequency'])
    plt.ylabel(parameter+' Error')
    if parameter in ['Frequency', 'FrequencyDerivative']:
        plt.ylabel('Relative '+parameter+' Error')
    plt.xlabel(x_parameter)
    plt.xscale('log')
    plt.yscale('log')
    # plt.legend(markerscale=4, loc = 'upper right')
    plt.savefig(SAVEPATH+'/Evaluation/SNRto'+parameter+'Error'+save_name+end_string)
    plt.show()

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

custom_lines = [plt.Line2D([0], [0], color=colors[0], lw=2, linestyle='dashed'),
                plt.Line2D([0], [0], color=colors[1], lw=2)]

correlation_list_flat = np.concatenate(correlation_list)
correlation_list_flat2 = np.concatenate(correlation_list2)


labels = {'EclipticLongitude': r'$\lambda$', 'EclipticLatitude': r'$\beta$','Frequency': '$f$ $($mHz$)$','FrequencyDerivative': r'$\dot{f}$ $ ($Hz/s$)$','Inclination': r'$\iota$','Amplitude': r'A', 'Polarization': r'$\psi$', 'InitialPhase': r'$\phi_0$'}

### histogram
# for parameter in parameters: 
# for parameter in  ['Amplitude']:
    # fig = plt.figure()
parameter_order = [
    "EclipticLatitude",
    "Frequency",
    "Amplitude",
    "InitialPhase",
    "EclipticLongitude",
    "FrequencyDerivative",
    "Inclination",
    "Polarization",
]

n_bins = 50
fig, axs = plt.subplots(2, 4, figsize=(10,6), constrained_layout=True)
for ax, parameter in zip(axs.flat, parameter_order):
    if parameter == 'Skylocation':
        ax.hist(np.abs(found_sources_matched_flat_df[parameter+'Error']), bins= np.logspace(-2,2, n_bins))
    elif parameter == 'Frequency':
        ax.hist(np.abs(found_sources_matched_flat_df[parameter+'Error']/pGB_injected_matched_flat_df[parameter]), bins= np.logspace(-8,-3.5, n_bins), log=False, density=False, histtype='step', label='ETH', linestyle='dashed')
        # ax.hist(np.abs(found_sources_matched_flat_df2[parameter+'Error']/pGB_injected_matched_flat_df2[parameter]), bins= np.logspace(-8,-3.5, n_bins), log=False, density=False, histtype='step', label='MM')
        # ax.hist(found_sources_matched_flat_df[parameter+'Error'], bins= np.linspace(-10**-7,10**-7, n_bins), log=True, density=True)
    elif parameter == 'FrequencyDerivative':
        ax.hist(np.abs(found_sources_matched_flat_df[parameter+'Error']), bins=np.logspace(-19,-13, n_bins), density=False, histtype='step', label='ETH', linestyle='dashed')
        # ax.hist(np.abs(found_sources_matched_flat_df2[parameter+'Error']), bins=np.logspace(-19,-13, n_bins), density=False, histtype='step', label='MM')
    elif parameter == 'Amplitude':
        ax.hist(found_sources_matched_flat_df[parameter+'Error'], bins=np.linspace(0,4,n_bins), density=False, histtype='step', label='ETH', linestyle='dashed')
        # ax.hist(found_sources_matched_flat_df2[parameter+'Error'], bins=np.linspace(0,4,n_bins), density=False, histtype='step', label='MM')
    else:
        ax.hist(found_sources_matched_flat_df[parameter+'Error'], bins=n_bins, density=False, histtype='step', label='ETH', linestyle='dashed')
        # ax.hist(found_sources_matched_flat_df2[parameter+'Error'], bins=n_bins, density=False, histtype='step', label='MM')
    ax.set_xlabel('$\Delta$'+parameter)

    if parameter == 'Skylocation':
        ax.set_xlabel('$\Delta$'+labels[parameter]+' (deg)')
        ax.set_xscale('log')
    if parameter == 'Amplitude':
        ax.set_xlabel(r'$\Delta \mathcal{A} / \mathcal{A}_{true}$')
        ax.set_xlim(0,4)
    if parameter == 'FrequencyDerivative':
        ax.set_xscale('log')
        ax.set_xlabel('$\Delta$'+labels[parameter])
        ax.set_xlim(10**-19,10**-13)
    if parameter == 'Frequency':
        ax.set_xscale('log')
        ax.set_xlim(10**-8,10**-3.5)
        # ax.ylim(0,10**3)
        ax.set_xlabel('$\Delta f / \mathcal{f}_{true}$')
    if parameter in ['EclipticLongitude', 'EclipticLatitude', 'Inclination', 'InitialPhase', 'Polarization']:
        ax.set_xlabel('$\Delta$'+labels[parameter]+' (rad)')
        ax.set_xlim(0,np.pi/2)
    ax.legend(custom_lines, ['ETH', 'MM'], loc='best')
    ax.grid(True)
    # ax.legend(loc='upper right')
# ax.yscale('log')
axs[0,1].legend(custom_lines, ['ETH', 'MM'], loc='upper left')
axs[1,1].legend(custom_lines, ['ETH', 'MM'], loc='upper left')
axs[0,0].set_ylabel('Count')
axs[1,0].set_ylabel('Count')
plt.savefig(SAVEPATH+'/Evaluation/_error_histogram'+save_name+end_string,dpi=300,bbox_inches='tight')
plt.show()

# customized settings
plot_parameter = {  # 'backend': 'ps',
    "font.family": "serif",
    "font.serif": "times",
    "font.size": 20,
    "axes.labelsize": "medium",
    "axes.titlesize": "medium",
    "legend.fontsize": "medium",
    "xtick.labelsize": "medium",
    "ytick.labelsize": "medium",
    "grid.color": "k",
    "grid.linestyle": ":",
    "grid.linewidth": 1,
    "savefig.dpi": 150,
}

# tell matplotlib about your param_plots
rcParams.update(plot_parameter)
# set nice figure sizes
fig_width_pt = 1.5*464.0  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (np.sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
ratio = golden_mean
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
fig_width = fig_width_pt * inches_per_pt  # width in inches
fig_height = fig_width * ratio  # height in inches
fig_size = [fig_width, fig_height]
fig_size_squared = [fig_width, fig_width]
rcParams.update({"figure.figsize": fig_size})


fig = plt.figure(figsize=fig_size, constrained_layout=True)
plt.hist(correlation_list_flat, bins= np.linspace(0,1,n_bins), histtype='step', label='ETH', linestyle='dashed', linewidth=2)
# plt.hist(correlation_list_flat2, bins= np.linspace(0,1,n_bins), histtype='step', label='MM', linewidth=2)
plt.xlabel(r'$\mathcal{M}$')
plt.ylabel('Count')
# plt.yscale('log')
plt.legend(custom_lines, ['ETH', 'MM'], loc='upper left')
plt.ylim(0,1100)
plt.xlim(0,1)
plt.grid(True)
plt.savefig(SAVEPATH+'/Evaluation/correlation_comparison'+save_name+end_string,dpi=300,bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=fig_size, constrained_layout=True)
plt.hist(found_sources_matched_flat_df['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', label='ETH', linestyle='dashed', linewidth=2, log=False)
# plt.hist(found_sources_matched_flat_df2['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', label='MM', linewidth=2, log=False)
plt.xlabel('Frequency')
plt.ylabel('Count')
plt.legend(custom_lines, ['ETH', 'MM'], loc='upper right')
plt.ylim(0,1000)
plt.xlim(0,0.03)
plt.grid(True)
plt.savefig(SAVEPATH+'/Evaluation/frequency_hist_comparison'+save_name+end_string,dpi=300,bbox_inches='tight')
plt.show()


found_sources_df = pd.concat([found_sources_matched_flat_df, found_sources_not_matched_flat_df])
found_sources_df2 = pd.concat([found_sources_matched_flat_df2, found_sources_not_matched_flat_df2])

fig = plt.figure(figsize=fig_size, constrained_layout=True)
counts, bins, bars = plt.hist(found_sources_df['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', label='ETH', linestyle='dashed', linewidth=2.5, log=False)
counts2, bins, bars = plt.hist(found_sources_df2['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', label='MM', linewidth=2, log=False)


fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True, figsize=fig_size )
# fig = plt.figure(figsize=fig_size, constrained_layout=True)
# counts, bins, bars = plt.hist(found_sources_df['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', label='ETH', linestyle='dashed', linewidth=2.5, log=False)
# counts2, bins, bars = plt.hist(found_sources_df2['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', label='MM', linewidth=2, log=False)
counts_matched, bins, bars = ax1.hist(found_sources_matched_flat_df['Frequency']*1000, bins= np.linspace(0,30,n_bins), histtype='step', label='ETH', linestyle='dashed', linewidth=2.5, log=False)
counts_matched2, bins, bars = ax1.hist(found_sources_matched_flat_df2['Frequency']*1000, bins= np.linspace(0,30,n_bins), histtype='step', label='MM', linewidth=2, log=False)
match_rate = counts_matched / counts
match_rate2 = counts_matched2 / counts2
# match_rate = np.nan_to_num(match_rate, nan=0.0)
# match_rate2 = np.nan_to_num(match_rate2, nan=0.0)
bins = np.linspace(0,30,n_bins)
center_bins = np.linspace((bins[1]+bins[0])/2,(bins[-2]+bins[-1])/2,n_bins-1)
ax2.set_xlabel('Frequency (mHz)')
ax1.set_ylabel('Count')
ax1.legend(custom_lines, ['ETH', 'MM'], loc='upper right')
ax1.set_ylim(0,1000)
ax1.set_xlim(0,30)
ax1.grid(True)
ax2.plot(center_bins,match_rate, '.', markersize=17, label='ETH')
ax2.plot(center_bins,match_rate2, 'P', markersize=6.3, label='MM')
ax2.set_ylabel('Match Rate')
ax2.grid(True)
ax2.legend(loc='center right')
plt.savefig(SAVEPATH+'/Evaluation/frequency_hist_all_found_comparison'+save_name+end_string,dpi=300,bbox_inches='tight')
plt.show()

print('end')