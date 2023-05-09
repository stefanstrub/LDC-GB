from re import A
# from matplotlib.lines import _LineStyle
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.font_manager
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



from Search import Search
from sources2 import *


# customized settings
plot_parameter = {  # 'backend': 'ps',
    "font.family" :'DeJavu Serif',
    "font.serif" : ["Computer Modern Serif"],
}


# customized settings
plot_parameter_big = {  # 'backend': 'ps',
    "font.family" :'DeJavu Serif',
    "font.serif" : ["Computer Modern Serif"],
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

SAVEPATH_sangria = grandparent+"/LDC/pictures/Sangria/"
save_name2 = 'Sangria_1year_dynamic_noise'
save_name0 = 'Sangria_12m'
save_name3 = 'Sangria_1_odd_dynamic_noise_SNR5'
# save_name3 = 'LDC1-4_2_optimized_second'
# save_name2 = 'Radler_1_full'
# save_name2 = 'Rxadler_1_full'
# save_name0 = 'LDC1-4_half_year'
save_name2 = 'Radler_12m'
save_name1 = 'Radler_24m'
# save_name = 'LDC1-4_half_year'
# save_name = 'Sangria_1_full_cut'

# duration = '3932160'
# duration = '7864320'
duration = '15728640'
# duration = '31457280'
save_name3 = 'Montana2022_'+duration
duration = '31457280'
# save_name1 = 'Montana2022_'+duration

save_names = [save_name0, save_name1, save_name2, save_name3]
SAVEPATHS = [SAVEPATH_sangria,SAVEPATH,SAVEPATH,SAVEPATH]

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

end_string = '_SNR_scaled_03_injected_snr5'
# end_string = 'correlation'
def load_files(save_path, save_name):
    found_sources_matched_flat_df = pd.read_pickle(save_path+'/found_sources_matched' +save_name+end_string+'_df')
    found_sources_not_matched_flat_df = pd.read_pickle(save_path+'/found_sources_not_matched' +save_name+end_string+'_df')
    pGB_injected_matched_flat_df = pd.read_pickle(save_path+'/injected_matched_windows' +save_name+end_string+'_df')
    pGB_injected_not_matched_flat_df = pd.read_pickle(save_path+'/injected_not_matched_windows' +save_name+end_string+'_df')
    match_list = np.load(save_path+'match_list' +save_name+end_string+'.npy', allow_pickle=True)
    pGB_best_list = np.load(save_path+'pGB_best_list' +save_name+end_string+'.npy', allow_pickle=True)
    match_best_list = np.load(save_path+'match_best_list' +save_name+end_string+'.npy', allow_pickle=True)

    pGB_best_list_flat = []
    for i in range(len(pGB_best_list)):
        for j in range(len(pGB_best_list[i])):
            pGB_best_list_flat.append(pGB_best_list[i][j])
    pGB_best_list_df = pd.DataFrame(np.asarray(pGB_best_list_flat))

    match_best_list_flat = np.concatenate(match_best_list)
    match_best_list_flat_array = {attribute: np.asarray([x[attribute] for x in match_best_list_flat]) for attribute in match_best_list_flat[0].keys()}
    match_best_list_flat_df = pd.DataFrame(match_best_list_flat_array)

    number_of_found_signals_not_matched = len(found_sources_not_matched_flat_df)
    number_of_matched_signals = len(found_sources_matched_flat_df)
    number_of_found_signals = number_of_matched_signals + number_of_found_signals_not_matched
    number_of_injected_signals = len(pGB_injected_matched_flat_df) + len(pGB_injected_not_matched_flat_df)

    print(number_of_matched_signals ,'matched signals out of', number_of_found_signals, 'found signals '+save_name)
    print('matched signals/found signals:', np.round(number_of_matched_signals/number_of_found_signals,2))
    print('number of injected signals:', np.round(number_of_injected_signals))
    return found_sources_matched_flat_df, found_sources_not_matched_flat_df, pGB_injected_matched_flat_df, pGB_injected_not_matched_flat_df, match_list, pGB_best_list_df, match_best_list_flat_df

found_sources_matched_list = []
found_sources_not_matched_list = []
pGB_injected_matched_list = []
pGB_injected_not_matched_list = []
match_list = []
pGB_best_list = []
match_best_list = []
for i, save_name in enumerate(save_names):
    found_sources_matched_flat_df, found_sources_not_matched_flat_df, pGB_injected_matched_flat_df, pGB_injected_not_matched_flat_df, match, pGB_best, match_best = load_files(SAVEPATHS[i], save_name)
    found_sources_matched_list.append(found_sources_matched_flat_df)
    found_sources_not_matched_list.append(found_sources_not_matched_flat_df)
    pGB_injected_matched_list.append(pGB_injected_matched_flat_df)
    pGB_injected_not_matched_list.append(pGB_injected_not_matched_flat_df)
    match_list.append(match)
    pGB_best_list.append(pGB_best)
    match_best_list.append(match_best)

# found_sources_matched_flat_df2, found_sources_not_matched_flat_df2, pGB_injected_matched_flat_df2, pGB_injected_not_matched_flat_df2, match_list2, pGB_best_list2, match_best_list2 = load_files(SAVEPATH, save_name2)
# found_sources_matched_flat_df3, found_sources_not_matched_flat_df3, pGB_injected_matched_flat_df3, pGB_injected_not_matched_flat_df3, match_list3, pGB_best_list3, match_best_list3 = load_files(SAVEPATH, save_name3)
# found_sources_matched_flat_df4, found_sources_not_matched_flat_df4, pGB_injected_matched_flat_df4, pGB_injected_not_matched_flat_df4, match_list4, pGB_best_list4, match_best_list4 = load_files(SAVEPATH, save_name4)

i = 0
SNR_threshold = 10
mask_matched = found_sources_matched_list[i]['IntrinsicSNR'] > SNR_threshold
mask_not_matched = found_sources_not_matched_list[i]['IntrinsicSNR'] > SNR_threshold
print(len(found_sources_matched_list[i][mask_matched]))
print(len(found_sources_not_matched_list[i][mask_not_matched]), len(found_sources_matched_list[i][mask_matched])/(len(found_sources_matched_list[i][mask_matched])+len(found_sources_not_matched_list[i][mask_not_matched])))
print(len(found_sources_matched_list[i]['Frequency'][found_sources_matched_list[i]['Frequency'] < 1.25*10**-3]))


#### plot SNR - frequency
markersize = 3
alpha = 0.5
i = 2
parameter_to_plot = 'IntrinsicSNR'
# parameter_to_plot = 'Amplitude'
fig = plt.figure()
# plt.plot(pGB_injectced_flat_df['Frequency']*10**3,pGB_injected_flat_df[parameter_to_plot], '.', color= colors[0], label = 'Injected', markersize= markersize, alpha = alpha)
# plt.plot(pGB_injected_matched_flat_df['Frequency']*10**3,pGB_injected_matched_flat_df[parameter_to_plot], '.', color= colors[1], label = 'Injected matched', markersize= markersize, alpha = alpha)
# plt.plot(pGB_injected_flat_df_high_SNR['Frequency']*10**3,pGB_injected_flat_df_high_SNR[parameter_to_plot],'.', color= colors[1], markersize= markersize, label = 'Injected SNR > 10', alpha = alpha)
plt.plot(found_sources_matched_list[i]['Frequency']*10**3,found_sources_matched_list[i][parameter_to_plot],'g.', label = 'Found matched', markersize= markersize, alpha = alpha, zorder = 5)
plt.plot(pGB_injected_not_matched_list[i]['Frequency']*10**3,pGB_injected_not_matched_list[i][parameter_to_plot], '+', color = 'r', label = 'Injected not matched', markersize= markersize, alpha = alpha)
plt.plot(found_sources_not_matched_list[i]['Frequency']*10**3,found_sources_not_matched_list[i][parameter_to_plot],'o',color= 'blue',  markerfacecolor='None', markersize= markersize, label = 'Found not matched', alpha = alpha)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('f (mHz)')
plt.xlim(0.3,100)
# plt.ylim(0.1,2000)
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

pGB_injected_list = [None] * 4
pGB_injected_high_SNR_list = [None] * 4
match_flat_list = [None] *4
for i in range(4):
    found_sources_matched_list[i] = get_errors(pGB_injected_matched_list[i], found_sources_matched_list[i])
    pGB_injected_list[i] = pd.concat([pGB_injected_matched_list[i], pGB_injected_not_matched_list[i]])
    pGB_injected_high_SNR_list[i] = pGB_injected_list[i][pGB_injected_list[i]['IntrinsicSNR'] > 10]
    match_flat_list[i] = np.concatenate(match_list[i])


##### error histogram

boundaries = {
    "EclipticLatitude": [-1.0, 1.0],
    "EclipticLongitude": [-np.pi, np.pi],
    "Frequency": [0.0003,0.033333],
    "Inclination": [-1.0, 1.0],
    "InitialPhase": [0.0, 2.0 * np.pi],
    "Polarization": [0.0, 1.0 * np.pi],
}

boundaries['GalactocentricLongitude'] = [0,360]
boundaries['GalactocentricLatitude'] = [-90,90]
# boundaries['Frequency'] = [frequencies_search[0][0],frequencies_search[-1][1]]
boundaries['IntrinsicSNR'] = [np.min(found_sources_matched_flat_df['IntrinsicSNR']),np.max(found_sources_matched_flat_df['IntrinsicSNR'])]
boundaries['FrequencyDerivative'] = [np.min(found_sources_matched_flat_df['FrequencyDerivative']),np.max(found_sources_matched_flat_df['FrequencyDerivative'])]
boundaries['Amplitude'] = [np.min(found_sources_matched_flat_df['Amplitude']),np.max(found_sources_matched_flat_df['Amplitude'])]
# parameter_x = 'GalacticLongitude'
# parameter_y = 'GalacticLatitude'
# parameter_x = 'EclipticLongitude'
# parameter_y = 'EclipticLatitude'
x_scale = 'log'
parameter_x = 'Frequency'
parameter_y = 'IntrinsicSNR'
n_bins = 50
x_coordinates = []
y_coordinates = []
for i in range(n_bins+1):
    length = (boundaries[parameter_x][1] - boundaries[parameter_x][0])/n_bins
    x_coordinates.append(boundaries[parameter_x][0]+length*i)
    length = (boundaries[parameter_y][1] - boundaries[parameter_y][0])/n_bins
    y_coordinates.append(boundaries[parameter_y][0]+length*i)

parameter_to_plot = 'SkylocationError'
parameter_to_plot = 'EclipticLatitudeError'
def get_2d_hist(sources, get_errors=True):
    error = []
    std = []
    count = []
    sources_parameter_x_sorted = sources.sort_values(by=parameter_x)
    # bin_boundaries_x = np.linspace(boundaries[parameter_x][0], boundaries[parameter_x][1],n_bins+1)
    # bin_boundaries_x = np.log10(bin_boundaries_x)
    bin_boundaries_y = np.linspace(boundaries[parameter_y][0], boundaries[parameter_y][1],n_bins+1)
    # bin_boundaries_y = np.log10(bin_boundaries_y)
    bin_boundaries_x = np.logspace(np.log10(boundaries[parameter_x][0]), np.log10(boundaries[parameter_x][1]),n_bins+1)
    bin_boundaries_y = np.logspace(np.log10(boundaries[parameter_y][0]), np.log10(boundaries[parameter_y][1]),n_bins+1)
    for i in range(n_bins):
        error.append([])
        std.append([])
        count.append([])
        # length = (boundaries[parameter_x][1] - boundaries[parameter_x][0])/n_bins
        start_index = np.searchsorted(sources_parameter_x_sorted[parameter_x],bin_boundaries_x[i], side='left')
        end_index = np.searchsorted(sources_parameter_x_sorted[parameter_x],bin_boundaries_x[i+1], side='left')
        section = sources_parameter_x_sorted[start_index:end_index]
        for j in range(n_bins):
            section_parameter_y_sorted = section.sort_values(by=parameter_y)
            # length = (boundaries[parameter_y][1] - boundaries[parameter_y][0])/n_bins
            start_index = np.searchsorted(section_parameter_y_sorted[parameter_y],bin_boundaries_y[j], side='left')
            end_index = np.searchsorted(section_parameter_y_sorted[parameter_y],bin_boundaries_y[j+1], side='left')
            field = section[start_index:end_index]
            change_to_deg = False
            if get_errors == True:
                if parameter_to_plot == 'SkylocationError':
                    change_to_deg = False
                if change_to_deg:
                    field[parameter_to_plot] *= 180/np.pi
                error[-1].append(np.mean(np.abs(field[parameter_to_plot])))
                std[-1].append(np.std(field[parameter_to_plot]))
            count[-1].append(len(field['Frequency']))
    for i in range(len(count)):
        for j in range(len(count[i])):
            if count[i][j] == 0:
                count[i][j] = np.nan
    return error, std, count, bin_boundaries_x, bin_boundaries_y
submission =0
error, std, count_not_matched, bin_boundaries_x, bin_boundaries_y = get_2d_hist(found_sources_not_matched_list[submission], get_errors=False)
error, std, count, bin_boundaries_x, bin_boundaries_y = get_2d_hist(found_sources_matched_list[submission], get_errors=True)

for i in range(len(count_not_matched)):
    for j in range(len(count_not_matched[i])):
        if np.isnan(count_not_matched[i][j]):
            count_not_matched[i][j] = 0

count_zeros = deepcopy(count)
for i in range(len(count_zeros)):
    for j in range(len(count_zeros[i])):
        if np.isnan(count_zeros[i][j]):
            count_zeros[i][j] = 0

match_ratio = np.asarray(count_zeros)/(np.asarray(count_zeros)+np.asarray(count_not_matched))

fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
fig.set_size_inches(8,10)
im = ax0.pcolormesh(bin_boundaries_x,bin_boundaries_y, np.array(error).T)
# im.set_clim(0,5)
fig.colorbar(im, ax=ax0)
ax0.set_title('mean '+parameter_to_plot)
im1 = ax1.pcolormesh(bin_boundaries_x,bin_boundaries_y, np.array(std).T)
# im1.set_clim(0,5)
fig.colorbar(im1, ax=ax1)
ax1.set_title('standard deviation')
im2 = ax2.pcolormesh(bin_boundaries_x,bin_boundaries_y, np.array(count).T)
fig.colorbar(im2, ax=ax2)
ax2.set_title('number of signals')
ax0.set_yscale('log')
ax1.set_yscale('log')
ax2.set_yscale('log')
ax0.set_xscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')
ax2.set_xlabel(parameter_x)
ax0.set_ylabel(parameter_y)
ax1.set_ylabel(parameter_y)
ax2.set_ylabel(parameter_y)
fig.tight_layout()
plt.savefig(SAVEPATH+'/Evaluation/'+parameter_x+parameter_y+parameter_to_plot+save_names[submission]+end_string)
plt.show()

fig, ax0 = plt.subplots(nrows=1, figsize=(10,6))
# fig.set_size_inches(6,6)
im = ax0.pcolormesh(bin_boundaries_x,bin_boundaries_y, np.array(match_ratio).T)
# im.set_clim(0,5)
cbar = fig.colorbar(im, ax=ax0)
# cbar.ax.set_ylabel('# of contacts', rotation=270, )
ax0.set_title('match ratio')
ax0.set_yscale('log')
ax0.set_xscale('log')
ax0.set_ylabel(parameter_y)
ax0.set_xlabel(parameter_x)
fig.tight_layout()
plt.savefig(SAVEPATH+'/Evaluation/'+parameter_x+parameter_y+'match_rate'+save_names[submission]+end_string)
plt.show()

##### plot SNR to error
# for parameter in  parameters +['Skylocation']:
# for parameter in  ['SkylocationError']:
x_parameter = 'Frequency'
# x_parameter = 'IntrinsicSNR'
i = 0
for parameter in  ['Frequency']:
    # x_parameter = parameter
    markersize = 4
    alpha = 0.5
    fig = plt.figure()
    # plt.plot(found_sources_matched_flat_df[parameter+'Error']/found_sources_matched_flat_df[parameter],found_sources_matched_flat_df[y_parameter],'.', label = 'found', markersize=markersize)
    # if parameter in ['Frequency', 'FrequencyDerivative']:
    #     plt.scatter(found_sources_matched_flat_df[x_parameter],found_sources_matched_flat_df[parameter+'Error']/found_sources_matched_flat_df['Frequency'], c=found_sources_matched_flat_df['Frequency'])
    # else:
    plt.plot(found_sources_matched_list[i][x_parameter],found_sources_matched_list[i][parameter+'Error'],'.', alpha =0.4, markersize=6)
        # plt.scatter(found_sources_matched_list[i][x_parameter],found_sources_matched_list[i][parameter], c=found_sources_matched_list[i]['Frequency'])
    plt.ylabel(parameter+' Error')
    if parameter in ['Frequency', 'FrequencyDerivative']:
        plt.ylabel('Relative '+parameter+' Error')
    plt.xlabel(x_parameter)
    plt.xscale('log')
    plt.yscale('log')
    # plt.legend(markerscale=4, loc = 'upper right')
    plt.savefig(SAVEPATH+'/Evaluation/SNRto'+parameter+'Error'+save_names[i]+end_string)
    plt.show()

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors[2] = '#00008B' 
colors[3] = '#BB5500' 

linestyle = ['solid', 'solid', 'dashed', 'solid']
# labels_plot = ['ETH 1 yr', 'MM 1 yr', 'ETH 0.5 yr', 'MM 0.5 yr']
# labels_plot = ['LDC1 0.5 yr', 'LDC1 1 yr', 'LDC2 1 yr', 'LDC1 2 yr']
labels_plot = ['LDC1 0.5 yr', 'LDC1 1 yr', 'LDC2 1 yr', 'LDC1 2 yr']
custom_lines = [plt.Line2D([0], [0], color=colors[0], lw=2, linestyle=linestyle[0]),
                plt.Line2D([0], [0], color=colors[1], lw=2, linestyle=linestyle[1]),
                plt.Line2D([0], [0], color=colors[2], lw=2, linestyle=linestyle[2]),
                plt.Line2D([0], [0], color=colors[3], lw=2, linestyle=linestyle[3])]


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

n_bins = 30
fig, axs = plt.subplots(2, 4, figsize=(10,6), constrained_layout=True)
for ax, parameter in zip(axs.flat, parameter_order):
    for i in range(4):
        if parameter == 'Skylocation':
            ax.hist(np.abs(found_sources_matched_flat_df[parameter+'Error']), bins= np.logspace(-2,2, n_bins))
        elif parameter == 'Frequency':
            ax.hist(np.abs(found_sources_matched_list[i][parameter+'Error']/pGB_injected_matched_list[i][parameter]), bins= np.logspace(-8,-3.5, n_bins), log=False, density=False, histtype='step', linestyle=linestyle[i], color=colors[i])
        elif parameter == 'FrequencyDerivative':
            ax.hist(np.abs(found_sources_matched_list[i][parameter+'Error']), bins=np.logspace(-19,-13, n_bins), density=False, histtype='step', linestyle=linestyle[i], color=colors[i])
        elif parameter == 'Amplitude':
            ax.hist(found_sources_matched_list[i][parameter+'Error'], bins=np.logspace(-4,1,n_bins), density=False, histtype='step', linestyle=linestyle[i], color=colors[i])
        elif parameter == 'Inclination':
            ax.hist(found_sources_matched_list[i][parameter+'Error'], bins=np.logspace(-5,np.log10(np.pi/2), n_bins), density=False, histtype='step', linestyle=linestyle[i], color=colors[i])
        elif parameter in ['EclipticLatitude', 'EclipticLongitude']:
            ax.hist(found_sources_matched_list[i][parameter+'Error'], bins=np.logspace(-5,np.log10(np.pi/2), n_bins), log= False, density=False, histtype='step', linestyle=linestyle[i], color=colors[i])
        else:
            ax.hist(found_sources_matched_list[i][parameter+'Error'], bins=n_bins, density=False, histtype='step', linestyle=linestyle[i], color=colors[i])
    ax.set_xlabel('$\Delta$'+parameter)
    if parameter == 'Skylocation':
        ax.set_xlabel('$\Delta$'+labels[parameter]+' (deg)')
        ax.set_xscale('log')
    if parameter == 'Inclination':
        ax.set_xlabel('$\Delta$'+labels[parameter]+' (rad)')
        ax.set_xlim(10**-5,np.pi/2)
        ax.set_xscale('log')
    if parameter == 'Amplitude':
        ax.set_xlabel(r'$\Delta \mathcal{A} / \mathcal{A}_{true}$')
        ax.set_xscale('log')
        ax.set_xlim(10**-4,10)
    if parameter == 'FrequencyDerivative':
        ax.set_xscale('log')
        ax.set_xlabel('$\Delta$'+labels[parameter])
        ax.set_xlim(10**-19,10**-13)
    if parameter == 'Frequency':
        ax.set_xscale('log')
        ax.set_xlim(10**-8,10**-3.5)
        # ax.ylim(0,10**3)
        ax.set_xlabel('$\Delta f / \mathcal{f}_{true}$')
    if parameter in [ 'InitialPhase', 'Polarization']:
        ax.set_xlabel('$\Delta$'+labels[parameter]+' (rad)')
        ax.set_xlim(0,np.pi/2)
    if parameter in ['EclipticLongitude', 'EclipticLatitude']:
        ax.set_xlabel('$\Delta$'+labels[parameter]+' (rad)')
        ax.set_xlim(10**-5,np.pi/2)
        ax.set_xscale('log')
    # ax.legend(custom_lines, [label1, label2, label3], loc='best')
    ax.grid(True)
    # ax.legend(loc='upper right')
# ax.yscale('log')
# axs[0,1].legend(custom_lines, [label1, label2, label3], loc='best')
# axs[1,1].legend(custom_lines, [label1, label2, label3], loc='best')
axs[0,0].legend(custom_lines, labels_plot, loc='upper left')
# axs[1,0].legend(custom_lines, [label1, label2, label3], loc='upper left')
# axs[0,2].legend(custom_lines, [label1, label2, label3], loc='upper left')
# axs[1,2].legend(custom_lines, [label1, label2, label3], loc='upper right')
axs[0,0].set_ylabel('Count')
axs[1,0].set_ylabel('Count')
plt.savefig(SAVEPATH+'/Evaluation/error_histogram'+save_names[0]+save_names[1]+save_names[2]+save_names[3]+end_string+'linear',dpi=300,bbox_inches='tight')
plt.show()


rcParams.update(plot_parameter_big)

upper_bound = 0.5
fig = plt.figure(figsize=fig_size, constrained_layout=True)
for i in range(4):
    plt.hist(match_flat_list[i], bins= np.linspace(0,upper_bound,n_bins), histtype='step', linestyle=linestyle[i], linewidth=2, color=colors[i])
plt.xlabel('$\delta$')
plt.ylabel('Count')
# plt.yscale('log')
plt.legend(custom_lines, labels_plot, loc='upper right')
# plt.ylim(0,2500)
plt.xlim(0,upper_bound)
plt.grid(True)
plt.savefig(SAVEPATH+'/Evaluation/correlation_comparison'+save_names[0]+save_names[1]+save_names[2]+save_names[3]+end_string,dpi=300,bbox_inches='tight')
plt.show()

upper_boundx = 0.5
upper_boundy = 50
n_bins = 100
vmax = [20,20,10,10]
for i in range(4):
    fig = plt.figure(figsize=fig_size, constrained_layout=True)
    plt.hist2d(match_flat_list[i], match_best_list[i]['IntrinsicSNR'], bins= (np.linspace(0,upper_boundx,n_bins),np.linspace(5,upper_boundy,n_bins)), vmin=0)
    plt.xlabel('$\delta$')
    plt.ylabel('SNR')
    plt.title(labels_plot[i])
    plt.colorbar( )
    plt.xlim(0,upper_boundx)
    plt.ylim(5,upper_boundy)
    plt.grid(True)
    plt.savefig(SAVEPATH+'/Evaluation/match_SNR_comparison'+save_names[i]+end_string,dpi=300,bbox_inches='tight')
    plt.show()



fig = plt.figure(figsize=fig_size, constrained_layout=True)
plt.hist(found_sources_matched_flat_df['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', linestyle=linestyle[i], linewidth=2, log=False)
# plt.hist(found_sources_matched_flat_df2['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', linewidth=2, log=False)
plt.xlabel('Frequency')
plt.ylabel('Count')
plt.legend(custom_lines,  labels_plot, loc='upper right')
plt.ylim(0,1000)
plt.xlim(0,0.03)
plt.grid(True)
plt.savefig(SAVEPATH+'/Evaluation/frequency_hist_comparison'+save_name+save_name2+end_string,dpi=300,bbox_inches='tight')
plt.show()

found_sources_df_list = [None]*4
for i in range(4):
    found_sources_df_list[i] = pd.concat([found_sources_matched_list[i], found_sources_not_matched_list[i]])
# found_sources_df2 = pd.concat([found_sources_matched_flat_df2, found_sources_not_matched_flat_df2])
# found_sources_df3 = pd.concat([found_sources_matched_flat_df3, found_sources_not_matched_flat_df3])
# found_sources_df4 = pd.concat([found_sources_matched_flat_df4, found_sources_not_matched_flat_df4])

# fig = plt.figure(figsize=fig_size)
# counts3, bins, bars = plt.hist(found_sources_df3['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', linewidth=2, log=False)
# counts, bins, bars = plt.hist(found_sources_df['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', linestyle=linestyle[i], linewidth=2.5, log=False)
# counts2, bins, bars = plt.hist(found_sources_df2['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', linewidth=2, log=False)


# fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True, figsize=fig_size )
# # fig = plt.figure(figsize=fig_size, constrained_layout=True)
# counts_matched, bins, bars = ax1.hist(found_sources_matched_flat_df['Frequency']*1000, bins= np.linspace(0,30,n_bins), histtype='step', linestyle=linestyle[i], linewidth=3, log=False)
# counts_matched2, bins, bars = ax1.hist(found_sources_matched_flat_df2['Frequency']*1000, bins= np.linspace(0,30,n_bins), histtype='step', linestyle=linestyle[i], linewidth=2.5, log=False)
# counts_matched3, bins, bars = ax1.hist(found_sources_matched_flat_df3['Frequency']*1000, bins= np.linspace(0,30,n_bins), histtype='step', linewidth=2, log=False)
# match_rate = counts_matched / counts
# match_rate2 = counts_matched2 / counts2
# match_rate3 = counts_matched3 / counts3
# # match_rate = np.nan_to_num(match_rate, nan=0.0)True)
# fig = plt.figure(figsize=fig_size)
# counts, bins, bars = plt.hist(found_sources_df['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=3, log=False, color=colors[0])
# counts2, bins, bars = plt.hist(found_sources_df2['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=2.5, log=False, color=colors[1])
# counts3, bins, bars = plt.hist(found_sources_df3['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[2])
# counts4, bins, bars = plt.hist(found_sources_df4['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[3])
# plt.xscale('log')
# plt.show()

n_bins = 30
match_rate_list = [None]*4
marker_list = ['.','o','P','P']
markersize_list = [20,7,8,7]
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, constrained_layout=True, figsize=[9,9])
for i in range(4):
    counts, bins, bars = ax1.hist(found_sources_df_list[i]['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=3, log=False, color=colors[i])
    # counts2, bins, bars = ax1.hist(found_sources_df2['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=2.5, log=False, color=colors[1])
    # counts3, bins, bars = ax1.hist(found_sources_df3['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[2])
    # counts4, bins, bars = ax1.hist(found_sources_df4['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[3])
    # fig = plt.figure(figsize=fig_size, constrained_layout=True)
    # counts, bins, bars = plt.hist(found_sources_df['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', linestyle=linestyle[i], linewidth=2.5, log=False)
    # counts2, bins, bars = plt.hist(found_sources_df2['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', linewidth=2, log=False)
    counts_matched, bins, bars = ax2.hist(found_sources_matched_list[i]['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=3, log=False, color=colors[i])
    # counts_matched2, bins, bars = ax2.hist(found_sources_matched_flat_df2['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=2.5, log=False, color=colors[1])
    # counts_matched3, bins, bars = ax2.hist(found_sources_matched_flat_df3['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[2])
    # counts_matched4, bins, bars = ax2.hist(found_sources_matched_flat_df4['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[3])
    match_rate_list[i] = counts_matched / counts
# match_rate2 = counts_matched2 / counts2
# match_rate3 = counts_matched3 / counts3
# match_rate4 = counts_matched4 / counts4
# match_rate = np.nan_to_num(match_rate, nan=0.0)
# match_rate2 = np.nan_to_num(match_rate2, nan=0.0)
bins = np.logspace(np.log10(0.3),np.log10(30),n_bins)
center_bins = np.logspace(np.log10((bins[1]+bins[0])/2),np.log10((bins[-2]+bins[-1])/2),n_bins-1)
ax3.set_xlabel('Frequency (mHz)')
ax1.set_ylabel('Submitted Signals')
ax2.set_ylabel('Matched Signals')
ax1.legend(custom_lines,  labels_plot, loc='upper right')
ax2.legend(custom_lines,  labels_plot, loc='upper right')
ax1.set_ylim(0,2500)
ax2.set_ylim(0,2500)
# ax3.set_ylim(0,1)
# ax1.set_xlim(0,30)
ax1.set_xscale('log')
ax2.set_xscale('log')
ax3.set_xscale('log')
ax1.grid(True)
ax2.grid(True)
for i in range(4):
    ax3.plot(center_bins, match_rate_list[i], marker_list[i],markersize= markersize_list[i],label=labels_plot[i], color=colors[i])
# ax3.plot(center_bins,match_rate3, '.', markersize=20, label=label3, color=colors[2])
# ax3.plot(center_bins,match_rate4, 'P', markersize=7, label=label4, color=colors[3])
# ax3.plot(center_bins,match_rate, 'o', markersize=8, label=label1 ,markerfacecolor='None', color=colors[0])
# ax3.plot(center_bins,match_rate2, 'P', markersize=7, label=label2, markerfacecolor='None', color=colors[1])
ax3.set_ylabel('Match Rate')
ax3.grid(True)
ax3.legend(loc='lower right')
plt.savefig(SAVEPATH+'/Evaluation/frequency_hist_all_found_comparison3'+save_names[0]+save_names[1]+save_names[2]+save_names[3]+end_string,dpi=300,bbox_inches='tight')
plt.show()

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, constrained_layout=True, figsize=[9,12])
for i in range(4):
    counts_injected, bins, bars = ax1.hist(pGB_injected_high_SNR_list[i]['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=3, log=False, color=colors[i])
    counts, bins, bars = ax2.hist(found_sources_df_list[i]['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=3, log=False, color=colors[i])
    # counts2, bins, bars = ax1.hist(found_sources_df2['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=2.5, log=False, color=colors[1])
    # counts3, bins, bars = ax1.hist(found_sources_df3['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[2])
    # counts4, bins, bars = ax1.hist(found_sources_df4['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[3])
    # fig = plt.figure(figsize=fig_size, constrained_layout=True)
    # counts, bins, bars = plt.hist(found_sources_df['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', linestyle=linestyle[i], linewidth=2.5, log=False)
    # counts2, bins, bars = plt.hist(found_sources_df2['Frequency'], bins= np.linspace(0,0.03,n_bins), histtype='step', linewidth=2, log=False)
    counts_matched, bins, bars = ax3.hist(found_sources_matched_list[i]['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=3, log=False, color=colors[i])
    # counts_matched2, bins, bars = ax2.hist(found_sources_matched_flat_df2['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linestyle=linestyle[i], linewidth=2.5, log=False, color=colors[1])
    # counts_matched3, bins, bars = ax2.hist(found_sources_matched_flat_df3['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[2])
    # counts_matched4, bins, bars = ax2.hist(found_sources_matched_flat_df4['Frequency']*1000, bins= np.logspace(np.log10(0.3),np.log10(30),n_bins), histtype='step', linewidth=2, log=False, color=colors[3])
    match_rate_list[i] = counts_matched / counts
# match_rate2 = counts_matched2 / counts2
# match_rate3 = counts_matched3 / counts3
# match_rate4 = counts_matched4 / counts4
# match_rate = np.nan_to_num(match_rate, nan=0.0)
# match_rate2 = np.nan_to_num(match_rate2, nan=0.0)
bins = np.logspace(np.log10(0.3),np.log10(30),n_bins)
center_bins = np.logspace(np.log10((bins[1]+bins[0])/2),np.log10((bins[-2]+bins[-1])/2),n_bins-1)
ax4.set_xlabel('Frequency (mHz)')
ax1.set_ylabel('Injected SNR > 10')
ax2.set_ylabel('Submitted')
ax3.set_ylabel('Matched')
ax1.legend(custom_lines,  labels_plot, loc='upper right')
# ax2.legend(custom_lines,  labels_plot, loc='upper right')
ax1.set_ylim(0,2500)
ax2.set_ylim(0,2500)
ax3.set_ylim(0,2500)
# ax3.set_ylim(0,1)
# ax1.set_xlim(0,30)
ax1.set_xscale('log')
ax2.set_xscale('log')
ax3.set_xscale('log')
ax4.set_xscale('log')
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
for i in range(4):
    ax4.plot(center_bins, match_rate_list[i], marker_list[i],markersize= markersize_list[i],label=labels_plot[i], color=colors[i])
# ax3.plot(center_bins,match_rate3, '.', markersize=20, label=label3, color=colors[2])
# ax3.plot(center_bins,match_rate4, 'P', markersize=7, label=label4, color=colors[3])
# ax3.plot(center_bins,match_rate, 'o', markersize=8, label=label1 ,markerfacecolor='None', color=colors[0])
# ax3.plot(center_bins,match_rate2, 'P', markersize=7, label=label2, markerfacecolor='None', color=colors[1])
ax4.set_ylabel('Match Rate')
ax4.grid(True)
ax4.legend(loc='lower right')
plt.savefig(SAVEPATH+'/Evaluation/frequency_hist_all_found_comparison4'+save_names[0]+save_names[1]+save_names[2]+save_names[3]+end_string,dpi=300,bbox_inches='tight')
plt.show()


def get_distance(amplitude, frequency, frequency_derivative):
    c = 3*10**8
    distance = 2/(96/5*np.pi**(8/3)*frequency**(11/3)/frequency_derivative)*np.pi*(2/3)*frequency**(2/3)/amplitude*c /3.086e+19 #to kpc
    return distance

def get_distance_for_dataframe(dataframe):
    dataframe['Distance'] = np.empty(len(dataframe['Amplitude']))
    postitive_fd_mask = dataframe['FrequencyDerivative'] >= 0
    distance = get_distance(dataframe['Amplitude'][postitive_fd_mask], dataframe['Frequency'][postitive_fd_mask],  dataframe['FrequencyDerivative'][postitive_fd_mask])
    dataframe['Distance'][postitive_fd_mask] = distance
    return dataframe

def get_galactic_coordinates(dataframe):
    coordinates_found = coord.SkyCoord(dataframe['EclipticLongitude'], dataframe['EclipticLatitude'], dataframe['Distance'], unit='rad', frame='barycentricmeanecliptic')
    dataframe['GalacticLongitude'] = coordinates_found.galactic.l.value
    dataframe['GalacticLatitude'] = coordinates_found.galactic.b.value
    dataframe['GalacticDistance'] = coordinates_found.galactic.distance.value
    return dataframe

# found_sources_matched_flat_df['Amplitude'] = pGB_injected_matched_flat_df['Amplitude']
# found_sources_matched_flat_df['FrequencyDerivative'] = pGB_injected_matched_flat_df['FrequencyDerivative']
for i in range(4):
    found_sources_matched_list[i] = get_distance_for_dataframe(found_sources_matched_list[i])
    found_sources_not_matched_list[i] = get_distance_for_dataframe(found_sources_not_matched_list[i])
    pGB_injected_matched_list[i] = get_distance_for_dataframe(pGB_injected_matched_list[i])
    # pGB_injected_not_matched_list[i] = get_distance_for_dataframe(pGB_injected_not_matched_list[i])

    found_sources_matched_list[i] = get_galactic_coordinates(found_sources_matched_list[i])
    found_sources_not_matched_list[i] = get_galactic_coordinates(found_sources_not_matched_list[i])
    pGB_injected_matched_list[i] = get_galactic_coordinates(pGB_injected_matched_list[i])
    # pGB_injected_not_matched_list[i] = get_galactic_coordinates(pGB_injected_not_matched_list[i])

i = 0
##### plot distance
markersize = 5
alpha = 0.5
fig = plt.figure()
postitive_fd_mask = found_sources_matched_list[i]['FrequencyDerivative'] >= 0
plt.plot(found_sources_matched_list[i]['GalacticLatitude'][postitive_fd_mask],found_sources_matched_list[i]['Distance'][postitive_fd_mask],'.', label = 'found', markersize=1)
# postitive_fd_mask = found_sources_not_matched_list[i]['FrequencyDerivative'] >= 0
# plt.plot(found_sources_not_matched_list[i]['GalacticLatitude'][postitive_fd_mask],found_sources_not_matched_list[i]['Distance'][postitive_fd_mask],'r.', label = 'Injected', markersize= 1, zorder=1)
postitive_fd_mask = pGB_injected_matched_list[i]['FrequencyDerivative'] >= 0
plt.plot(pGB_injected_matched_list[i]['GalacticLatitude'][postitive_fd_mask],pGB_injected_matched_list[i]['Distance'][postitive_fd_mask],'g.', label = 'Injected', markersize= 1, zorder=1)
# postitive_fd_mask = pGB_injected_not_matched_list[i]['FrequencyDerivative'] >= 0
# plt.plot(pGB_injected_not_matched_list[i]['GalacticLatitude'][postitive_fd_mask],pGB_injected_not_matched_list[i]['Distance'][postitive_fd_mask],'+', label = 'not matched', color = 'r', markersize=2, zorder= 1)
# postitive_fd_mask = pGB_injected_flat_highSNR_df['FrequencyDerivative'] >= 0
# plt.plot(pGB_injected_flat_highSNR_df['GalacticLatitude'][postitive_fd_mask],pGB_injected_flat_highSNR_df['Distance'][postitive_fd_mask],'+', label = 'injected SNR>10', color = 'r', markersize=2, zorder= 4)
plt.xlabel('GalacticLatitude')
plt.ylabel('Distance [kpc]')
plt.legend(markerscale=4, loc = 'upper right')
plt.savefig(SAVEPATH+'/galacticLatitudeDistance'+save_names[i])
plt.show()



print('end')