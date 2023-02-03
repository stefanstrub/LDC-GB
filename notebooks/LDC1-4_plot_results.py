#%%
from re import A
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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
import yaml

from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries, window
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr


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
fig_width_pt = 1.5*464.0  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (np.sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
ratio = golden_mean
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
fig_width = fig_width_pt * inches_per_pt  # width in inches
fig_height = fig_width * ratio  # height in inches
fig_size = [fig_width, fig_height]
fig_size_squared = [fig_width, fig_width]
rcParams.update({"figure.figsize": fig_size})

# get current directory
path = os.getcwd()
 
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)

DATAPATH = "/home/stefan/LDC/Radler/data"
DATAPATH = grandparent+"/LDC/Radler/data"
SAVEPATH = grandparent+"/LDC/Radler/LDC1-4_evaluation"

data_file_name = 'Montana'
# data_file_name = 'APC'
# data_file_name = 'ETH'
save_name = 'LDC1-4_4mHz'+data_file_name

found_sources_matched_flat_df = {}
found_sources_not_matched_flat_df = {}
pGB_injected_not_matched_flat_df = {}
pGB_injected_matched_flat_df = {}
for data_file_name in ['Montana', 'ETH', 'APC',]:
    save_name = 'LDC1-4_4mHz'+data_file_name
    found_sources_matched_flat_df[data_file_name] = pd.read_pickle(SAVEPATH+'/found_sources_matched' +save_name+'.pkl')  
    found_sources_not_matched_flat_df[data_file_name] = pd.read_pickle(SAVEPATH+'/found_sources_not_matched' +save_name+'.pkl')  
    pGB_injected_not_matched_flat_df[data_file_name] = pd.read_pickle(SAVEPATH+'/injected_not_matched' +save_name+'.pkl')  
    pGB_injected_matched_flat_df[data_file_name] = pd.read_pickle(SAVEPATH+'/injected_matched' +save_name+'.pkl')  


#### plot SNR - frequency
markersize = 5
alpha = 0.4
parameter_to_plot = 'IntrinsicSNR'
fig = plt.figure()
for data_file_name, color in [('Montana','green'), ('ETH', 'blue'), ('APC', 'red')]:
    print(data_file_name, color)
    # plt.plot(pGB_injected_flat_df[data_file_name]['Frequency']*10**3,pGB_injected_flat_df[data_file_name]['IntrinsicSNR'], '.', color= color, label = 'Injected', markersize= markersize, alpha = alpha)
    # plt.plot(pGB_injected_matched_flat_df[data_file_name]['Frequency']*10**3,pGB_injected_matched_flat_df[data_file_name]['IntrinsicSNR'], '.', color= color, label = 'Injected matched', markersize= markersize, alpha = alpha)
    # plt.plot(pGB_injected_flat_df_high_SNR[data_file_name]['Frequency']*10**3,pGB_injected_flat_df_high_SNR[data_file_name]['IntrinsicSNR'],'.', color= color, markersize= markersize, label = 'Injected SNR > 10', alpha = alpha)
    plt.plot(found_sources_matched_flat_df[data_file_name]['Frequency']*10**3,found_sources_matched_flat_df[data_file_name]['IntrinsicSNR'],'.', color= color, label = 'Found matched ' + data_file_name, markersize= markersize, alpha = alpha)
    plt.plot(pGB_injected_not_matched_flat_df[data_file_name]['Frequency']*10**3,pGB_injected_not_matched_flat_df[data_file_name]['IntrinsicSNR'], '+', color= color, label = 'Injected not matched ' + data_file_name, markersize= markersize, alpha = alpha)
    plt.plot(found_sources_not_matched_flat_df[data_file_name]['Frequency']*10**3,found_sources_not_matched_flat_df[data_file_name]['IntrinsicSNR'],'o',color= color,  markerfacecolor='None', markersize= markersize, label = 'Found not matched ' + data_file_name, alpha = alpha)
plt.yscale('log')
# plt.xscale('log')
plt.xlabel('f (mHz)')
# plt.xlim(0.3,30)
# plt.ylim(0.1,2000)
if parameter_to_plot == 'IntrinsicSNR':
    plt.ylabel(r'$\mathrm{SNR}_{\mathrm{opt}}$')
else:
    plt.ylabel(parameter_to_plot)    
plt.legend(markerscale=3, loc = 'upper right')
# plt.savefig(SAVEPATH+'/Evaluation/'+parameter_to_plot+save_name+'injected_all_found',dpi=300,bbox_inches='tight')
plt.show()


#### plot SNR - frequency
markersize = 7
alpha = 0.4
parameter_to_plot = 'IntrinsicSNR'
fig = plt.figure()
for data_file_name, color in [('Montana','green'), ('ETH', 'blue'), ('APC', 'red')]:
    print(data_file_name, color)
    # plt.plot(pGB_injected_flat_df[data_file_name]['Frequency']*10**3,pGB_injected_flat_df[data_file_name]['IntrinsicSNR'], '.', color= color, label = 'Injected', markersize= markersize, alpha = alpha)
    # plt.plot(pGB_injected_matched_flat_df[data_file_name]['Frequency']*10**3,pGB_injected_matched_flat_df[data_file_name]['IntrinsicSNR'], '.', color= color, label = 'Injected matched', markersize= markersize, alpha = alpha)
    # plt.plot(pGB_injected_flat_df_high_SNR[data_file_name]['Frequency']*10**3,pGB_injected_flat_df_high_SNR[data_file_name]['IntrinsicSNR'],'.', color= color, markersize= markersize, label = 'Injected SNR > 10', alpha = alpha)
    plt.plot(found_sources_matched_flat_df[data_file_name]['Frequency']*10**3,found_sources_matched_flat_df[data_file_name]['IntrinsicSNR'],'.', color= color, label = 'Found matched ' + data_file_name, markersize= markersize, alpha = alpha)
    plt.plot(pGB_injected_not_matched_flat_df[data_file_name]['Frequency']*10**3,pGB_injected_not_matched_flat_df[data_file_name]['IntrinsicSNR'], '+', color= color, label = 'Injected not matched ' + data_file_name, markersize= markersize, alpha = alpha)
    plt.plot(found_sources_not_matched_flat_df[data_file_name]['Frequency']*10**3,found_sources_not_matched_flat_df[data_file_name]['IntrinsicSNR'],'o',color= color,  markerfacecolor='None', markersize= markersize, label = 'Found not matched ' + data_file_name, alpha = alpha)
plt.yscale('log')
# plt.xscale('log')
plt.xlabel('f (mHz)')
# plt.xlim(0.3,30)
plt.ylim(1,1000)
if parameter_to_plot == 'IntrinsicSNR':
    plt.ylabel(r'$\mathrm{SNR}_{\mathrm{opt}}$')
else:
    plt.ylabel(parameter_to_plot)  
legend_elements = [mlines.Line2D([], [], color='grey', marker='.', linestyle='None',
                          markersize=10, label='Found matched'),
                   mlines.Line2D([], [], color='grey', marker='+', linestyle='None',
                          markersize=10, label='Injected not matched'),
                   mlines.Line2D([], [], color='grey', marker='o', markerfacecolor='None', linestyle='None',
                          markersize=10, label='Found not matched')]
list_lab = ['Found matched','Injected not matched','Found not matched']
plt.legend(legend_elements,list_lab, loc='lower right')
plt.savefig(SAVEPATH+'/Evaluation/'+parameter_to_plot+save_name+'injected_all_found',dpi=300,bbox_inches='tight')
plt.show()

print('end')