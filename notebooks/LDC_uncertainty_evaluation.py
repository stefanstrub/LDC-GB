from re import A
# from matplotlib.lines import _LineStyle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LogNorm
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
import itertools
from KDEpy import FFTKDE
import pickle

from scipy import integrate



# customized settings
plot_parameter = {  # 'backend': 'ps',
    "font.family" :'DeJavu Serif',
    "font.serif" : ["Computer Modern Serif"],
    "mathtext.fontset": "cm",
}


# customized settings
plot_parameter_big = {  # 'backend': 'ps',
    "font.family" :'DeJavu Serif',
    "font.serif": "Times",
    # "font.serif" : ["Computer Modern Serif"],
    "font.size": 20,
    "mathtext.fontset": "cm",
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
rcParams.update(plot_parameter_big)
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
SAVEPATH = grandparent+"/LDC/pictures/LDC1-4/"
SAVEPATH = grandparent+"/LDC/Radler/LDC1-4_evaluation/"
SAVEPATH_sangria = grandparent+"/LDC/Sangria/evaluation/"
SAVEPATH_sangria = grandparent+"/LDC/pictures/Sangria/"

save_name2 = 'Sangria_12m'
save_name4 = 'Radler_6m'
save_name3 = 'Radler_12m'
save_name1 = 'Radler_24m'

# duration = '3932160'
# duration = '7864320'
duration = '15728640'
# duration = '31457280'
# save_name4 = 'Montana2022_'+duration
duration = '31457280'


save_names = [save_name1, save_name2, save_name3, save_name4]
SAVEPATHS = [SAVEPATH,SAVEPATH_sangria,SAVEPATH,SAVEPATH]

save_names = [save_name2]
SAVEPATHS = [SAVEPATH_sangria]

# data_set = 1
chain_paths = []
for data_set in range(len(save_names)):
    chain_paths.append(SAVEPATHS[data_set] + 'Chains_gpu_partial07_'+save_names[data_set]+'/')


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

parameters = [
    "Amplitude",
    "EclipticLatitude",
    "EclipticLongitude",
    "Frequency",
    "FrequencyDerivative",
    "Inclination"
]

labels = {'EclipticLongitude': r'$\lambda$'+' (rad)', 'EclipticLatitude': r'$\beta$'+' (rad)','Frequency': r'$f / f_\mathrm{true}$','FrequencyDerivative': r'$\dot{f}$ $ ($Hz/s$)$','Inclination': r'$\iota$'+' (rad)','Amplitude': r'$ \mathcal{A}$', 'Polarization': r'$\psi$'+' (rad)', 'InitialPhase': r'$\phi_0$'+' (rad)'}

end_string = '_SNR_scaled_03_injected_snr5'
# end_string = 'correlation'
def load_files(save_path, save_name):
    found_sources_matched_flat_df = pd.read_pickle(save_path+'/found_sources_matched_' +save_name+end_string+'_df')
    found_sources_not_matched_flat_df = pd.read_pickle(save_path+'found_sources_not_matched_' +save_name+end_string+'_df')
    pGB_injected_matched_flat_df = pd.read_pickle(save_path+'/injected_matched_windows_' +save_name+end_string+'_df')
    pGB_injected_not_matched_flat_df = pd.read_pickle(save_path+'/injected_not_matched_windows_' +save_name+end_string+'_df')
    match_list = np.load(save_path+'match_list_' +save_name+end_string+'.npy', allow_pickle=True)
    pGB_best_list = np.load(save_path+'pGB_best_list_' +save_name+end_string+'.npy', allow_pickle=True)
    match_best_list = np.load(save_path+'match_best_list_' +save_name+end_string+'.npy', allow_pickle=True)

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

def normalize(array, min, max):
    return (array-min)/(max-min)

def check_if_in_confidence_region_1D_std(parameter_x):
    std = np.std(df[parameter_x].values)
    mean = np.mean(df[parameter_x].values)
    true_parameter = pGB_injected_matched_list[data_set][parameter_x][index_of_signal]
    if true_parameter > mean-std and true_parameter < mean+std:
        in_contour = True
    else:
        in_contour = False
    return in_contour

def check_if_in_confidence_region_1D_kde(parameter_x):
    min_x = np.min(df[parameter_x])
    max_x = np.max(df[parameter_x])
    # Notice how bw (standard deviation), kernel, weights and grid points are set
    x, pdf = FFTKDE(bw=0.1, kernel='gaussian').fit(np.array(normalize(df[parameter_x], min_x, max_x)), weights=None).evaluate(2**8)

    hist_indexes_sorted = np.unravel_index(np.argsort(-pdf, axis=None), pdf.shape)
    pdf_sorted = pdf[hist_indexes_sorted]*(x[1]-x[0])

    # plt.figure()
    # plt.plot(x, pdf); plt.tight_layout()
    # plt.show()

    # plt.figure()
    # plt.hist(df[parameter_x].values)
    # plt.show()

    sum = 0
    for i in range(len(pdf_sorted)):
        sum += pdf_sorted[i]
        if sum > confidence_threshold:
            largest_index_in_contour = i
            break

    contour = np.zeros_like(pdf)
    for i in range(largest_index_in_contour):
        contour[hist_indexes_sorted[0][i]] = 1
    true_parameter = pGB_injected_matched_list[data_set][parameter_x][index_of_signal]
    index_true = np.searchsorted(x, normalize(true_parameter, min_x, max_x))

    if index_true > len(pdf)-1:
        index_true -= 1

    if index_true == 0:
        in_contour = False
        return in_contour

    if index_true > len(pdf)-1:
        in_contour = False
        return in_contour
    
    if contour[index_true] == 1:
        in_contour = True
    else:
        in_contour = False

    return in_contour

def check_if_in_confidence_region_1D(parameter_x):
    number_of_samples = len(df[parameter_x].values)
    hist_both = np.histogram(df[parameter_x].values, bins=50)
    # plt.figure()
    # plt.plot(hist_both[1][:-1], hist_both[0])
    # plt.show()

    hist = hist_both[0]/number_of_samples
    hist_axis = hist_both[1]
    hist_indexes_sorted = np.unravel_index(np.argsort(-hist, axis=None), hist.shape)
    hist_sorted = hist[hist_indexes_sorted]

    sum = 0
    for i in range(len(hist_sorted)):
        sum += hist_sorted[i]
        if sum > confidence_threshold:
            largest_index_in_contour = i
            break

    contour = np.zeros_like(hist)
    for i in range(largest_index_in_contour):
        contour[hist_indexes_sorted[0][i]] = 1
    true_parameter = pGB_injected_matched_list[data_set][parameter_x][index_of_signal]
    index_true = np.searchsorted(hist_axis, true_parameter)

    if index_true > len(hist)-1:
        index_true -= 1

    if index_true == 0:
        in_contour = False
        return in_contour

    if index_true > len(hist)-1:
        in_contour = False
        return in_contour
    
    if contour[index_true] == 1:
        in_contour = True
    else:
        in_contour = False

    return in_contour

def check_if_in_confidence_region_2D(parameter_x, parameter_y):
    number_of_samples = len(df[parameter_x].values)
    hist_values_axes = np.histogram2d(df[parameter_x].values, df[parameter_y].values, bins=10)
    # plt.figure()
    # plt.plot(hist_values_axes[1][:-1], hist_values_axes[0])
    # plt.show()

    hist = hist_values_axes[0]/number_of_samples
    hist_indexes_sorted = np.unravel_index(np.argsort(-hist, axis=None), hist.shape)
    hist_sorted = hist[hist_indexes_sorted]

    sum = 0
    for i in range(len(hist_sorted)):
        sum += hist_sorted[i]
        if sum > confidence_threshold:
            largest_index_in_contour = i
            break

    contour = np.zeros_like(hist)
    for i in range(largest_index_in_contour):
        contour[hist_indexes_sorted[0][i],hist_indexes_sorted[1][i]] = 1
    true_parameter_x = pGB_injected_matched_list[data_set][parameter_x][index_of_signal]
    true_parameter_y = pGB_injected_matched_list[data_set][parameter_y][index_of_signal]
    index_true = [np.searchsorted(hist_values_axes[1], true_parameter_x), np.searchsorted(hist_values_axes[data_set], true_parameter_y)]

    if index_true[0] > len(hist)-1:
        index_true[0] -= 1
    if index_true[1] > len(hist)-1:
        index_true[1] -= 1


    if index_true[0] > len(hist)-1 or index_true[1] > len(hist)-1:
        in_contour = False
        return in_contour
    
    if index_true[0] == 0 or index_true[1] == 0:
        in_contour = False
        return in_contour

    if contour[index_true[1],index_true[0]] == 1:
        in_contour = True
    else:
        in_contour = False

    # plt.figure()
    # plt.imshow(contour)
    # plt.xlabel(parameter_x)
    # plt.ylabel(parameter_y)
    # plt.show()

    # plt.figure()
    # plt.hist2d(df[parameter_x].values, df[parameter_y].values, bins=20)
    # plt.xlabel(parameter_x)
    # plt.ylabel(parameter_y)
    # plt.show()


    return in_contour

def check_if_in_confidence_region_2d_kde(parameter_x, parameter_y):
    min_x = np.min(df[parameter_x])
    max_x = np.max(df[parameter_x])
    min_y = np.min(df[parameter_y])
    max_y = np.max(df[parameter_y])

    grid_points = 2**6

    # grid_points = 2**5+1
    # ax = np.linspace(-0.15,1.15,grid_points)
    # mypdf,axes = fastKDE.pdf(normalize(df[parameter_x], min_x, max_x),normalize(df[parameter_y], min_y, max_y), axes=[ax,ax])
    # x = ax

    data = np.array([normalize(df[parameter_x], min_x, max_x),normalize(df[parameter_y], min_y, max_y)]).T
    grid, mypdf = FFTKDE(bw=0.1, kernel='gaussian').fit(data).evaluate(grid_points)
    x, y = np.unique(grid[:, 0]), np.unique(grid[:, 1])
    mypdf = mypdf.reshape(grid_points, grid_points).T


    # Plot the kernel density estimate
    # N = 16
    # plt.figure()
    # plt.contour(x, y, z, N, linewidths=0.8, colors='k')
    # plt.contourf(x, y, z, N, cmap="RdBu_r")
    # plt.plot(data[:, 0], data[:, 1], 'ok', ms=3)
    # plt.show()

    # while np.min(mypdf) < 0:
    #     ax_expand = (np.random.rand(1)*0.2+0.1)[0]
    #     ax = np.linspace(-ax_expand,1+ax_expand,numPoints)
    #     mypdf,axes = fastKDE.pdf(normalize(df[parameter_x], min_x, max_x),normalize(df[parameter_y], min_y, max_y), axes=[ax,ax])

    # plt.figure()
    # plt.imshow(mypdf)
    # plt.xlabel(parameter_x)
    # plt.ylabel(parameter_y)
    # plt.show()

    # sum = np.sum(mypdf)*(ax[1]-ax[0])**2

    mypdf_indexes_sorted = np.unravel_index(np.argsort(-mypdf, axis=None), mypdf.shape)
    mypdf_sorted = mypdf[mypdf_indexes_sorted]*(x[1]-x[0])**2

    sum = 0
    for i in range(len(mypdf_sorted)):
        sum += mypdf_sorted[i]
        if sum > confidence_threshold:
            largest_index_in_contour = i
            break

    contour = np.zeros_like(mypdf)
    for i in range(largest_index_in_contour):
        contour[mypdf_indexes_sorted[0][i],mypdf_indexes_sorted[1][i]] = 1
    x_normalized = normalize(pGB_injected_matched_list[data_set][parameter_x][index_of_signal], min_x, max_x)
    y_normalized = normalize(pGB_injected_matched_list[data_set][parameter_y][index_of_signal], min_y, max_y)
    index_true = [np.searchsorted(x, x_normalized), np.searchsorted(x, y_normalized)]
    if index_true[0] > grid_points-1:
        index_true[0] -= 1
    if index_true[1] > grid_points-1:
        index_true[1] -= 1

    if contour[index_true[1],index_true[0]] == 1:
        in_contour = True
    else:
        in_contour = False
    # print('In contour ',true_in_contour)
    # contour[index_true[1],index_true[0]] = 0.7

    # plt.figure()
    # plt.imshow(contour)
    # plt.xlabel(parameter_x)
    # plt.ylabel(parameter_y)
    # plt.show()

    return in_contour


# frequencies_in_folder = []
# for filename in onlyfiles:
#     frequencies_in_folder.append(float(filename[9:filename.find('pHz')])/10**12)
# frequencies_in_folder = np.array(frequencies_in_folder)
# frequencies_in_folder = np.sort(frequencies_in_folder)


in_contour_f_fd_list = []
in_contour_f_ip_list = []
in_contour_sky_list = []
in_contour_amp_inclination_list = []
in_contour_list = []
in_contour = {}
std = {}
mean = {}
for parameter_x, parameter_y in itertools.combinations(parameters, 2):
    in_contour[parameter_x,parameter_y] = []

for parameter_x in parameters:
    in_contour[parameter_x,parameter_x] = []
for parameter in parameters:
    in_contour[parameter] = []
    mean[parameter] = []
    std[parameter] = []

confidence_threshold = 0.68

# samples = np.random.normal(0.8, 0.01, 1000)
# std = np.std(samples)
# mean = np.mean(samples)
# print(mean, std)
# print(np.arcsin(mean), np.arcsin(std))
# samples_arc = np.arcsin(samples)
# std_arc = np.std(samples_arc)
# mean_arc = np.mean(samples_arc)
# print(mean_arc, std_arc)




# print(np.sin(mean), np.sin(std))

# for index, filename in enumerate(onlyfiles):
#     if index % int(len(onlyfiles)/10) == 0:
#         print(np.round(index/len(onlyfiles),2))
#     df = pd.read_csv(chain_path+filename)
#     for parameter in parameters:
#         mean[parameter].append(np.mean(df[parameter]))
#         std[parameter].append(np.std(df[parameter]))

# found_sources_list = []
# mean_list = []
# std_list = []

# for data_set in range(len(save_names)):
#     onlyfiles = [f for f in os.listdir(chain_paths[data_set]) if os.path.isfile(os.path.join(chain_paths[data_set], f))]
#     found_sources = pd.concat([found_sources_matched_list[data_set], found_sources_not_matched_list[data_set]])
#     found_sources = found_sources.sort_values(by=['Frequency'])
#     found_sources = found_sources.reset_index()
#     std = {}
#     mean = {}
#     for parameter in parameters:
#         mean[parameter] = []
#         std[parameter] = []
#     for index in range(len(onlyfiles)):
#         # if index > 46:
#         #     continue
#         # index_of_signal = np.searchsorted(found_sources['Frequency'], frequencies_in_folder[index])
#         # if index_of_signal == len(onlyfiles): 
#         #     index_of_signal = index_of_signal-1
#         # if index_of_signal < len(onlyfiles)-1: 
#         #     if abs(found_sources['Frequency'][index_of_signal]-frequencies_in_folder[index]) > abs(found_sources['Frequency'][index_of_signal+1]-frequencies_in_folder[index]):
#         #         index_of_signal = index_of_signal+1
#         # if index_of_signal != 0: 
#         #     if abs(found_sources['Frequency'][index_of_signal]-frequencies_in_folder[index]) > abs(found_sources['Frequency'][index_of_signal-1]-frequencies_in_folder[index]):
#         #         index_of_signal = index_of_signal-1
#         index_of_signal = index
#         # index_of_signal = 4804
#         # index_of_signal = 502
#         # index_of_signal = 9
#         if index % int(len(onlyfiles)/10) == 0:
#             print(np.round(index/len(onlyfiles),2))
#         # print('frequency'+str(int(np.round(found_sources['Frequency'][index_of_signal]*10**9)))+'Sangria_1_full_cut.csv')
#         try:
#             df = pd.read_csv(chain_paths[data_set]+'frequency'+str(int(np.round(found_sources['Frequency'][index_of_signal]*10**12)))+'pHz'+ save_names[data_set] +'.csv')
#         except:
#             pass
        
#         df['Frequency'] *= 10**-3
#         df['EclipticLatitude'] = np.sin(df['EclipticLatitude'])
#         for parameter in parameters:
#             mean[parameter].append(np.mean(df[parameter]))
#             std[parameter].append(np.std(df[parameter]))
#     found_sources_list.append(found_sources)
#     mean_list.append(mean)
#     std_list.append(std)


# pickle.dump(mean_list, open(SAVEPATH+'/mean_list' +save_name+end_string+'.npy', "wb"))
# pickle.dump(std_list, open(SAVEPATH+'/std_list' +save_name+end_string+'.npy', "wb"))
# pickle.dump(found_sources_list, open(SAVEPATH+'/found_sources_list' +save_name+end_string+'.npy', "wb"))

mean_list = pickle.load(open(SAVEPATH+'/mean_list' +save_name+end_string+'.npy', "rb"))
std_list = pickle.load(open(SAVEPATH+'/std_list' +save_name+end_string+'.npy', "rb"))
found_sources_list = pickle.load(open(SAVEPATH+'/found_sources_list' +save_name+end_string+'.npy', "rb"))


for data_set in range(len(save_names)):
    for i in range(len(mean_list[data_set]['EclipticLongitude'])):
        if mean_list[data_set]['EclipticLongitude'][i] < 0:
            mean_list[data_set]['EclipticLongitude'][i] += 2*np.pi


def angular_distance(ra1, dec1, ra2, dec2):
    return np.arccos(np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2))


def f(x, y):
    return np.abs(np.sin(x))


prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
blue = deepcopy(colors[0])
colors[0] = deepcopy(colors[1])
colors[1] = deepcopy(blue)
colors[2] = '#00008B' 
# colors[3] = '#BB5500' 

linestyle = ['solid', 'dashed', 'solid', 'solid']
# labels_plot = ['ETH 1 yr', 'MM 1 yr', 'ETH 0.5 yr', 'MM 0.5 yr']
# labels_plot = ['LDC1 0.5 yr', 'LDC1 1 yr', 'LDC2 1 yr', 'LDC1 2 yr']
# labels_plot = ['LDC1 2 yr', 'LDC2 1 yr', 'LDC1 1 yr', 'LDC1 0.5 yr']
labels_plot = ['Radler 2 yr', 'Sangria 1 yr', 'Radler 1 yr', 'Radler 0.5 yr']
line_width = 3
custom_lines = [plt.Line2D([0], [0], color=colors[0], lw=line_width, linestyle=linestyle[0]),
                plt.Line2D([0], [0], color=colors[1], lw=line_width, linestyle=linestyle[1]),
                plt.Line2D([0], [0], color=colors[2], lw=line_width, linestyle=linestyle[2]),
                plt.Line2D([0], [0], color=colors[3], lw=line_width, linestyle=linestyle[3])]


line_width = 3
data_set = 2
target_frequency = 0.004169906
index = np.searchsorted(found_sources_list[data_set]['Frequency'], target_frequency)
found_sources_list[data_set]['Frequency'][index]
if abs(found_sources_list[data_set]['Frequency'][index]-target_frequency) > abs(found_sources_list[data_set]['Frequency'][index-1]-target_frequency):
    index = index-1
    
print(mean_list[data_set]['Frequency'][index], mean_list[data_set]['EclipticLatitude'][index], mean_list[data_set]['EclipticLongitude'][index])
print(mean_list[data_set]['EclipticLatitude'][index]-std_list[data_set]['EclipticLatitude'][index],mean_list[data_set]['EclipticLatitude'][index]+std_list[data_set]['EclipticLatitude'][index])
print(mean_list[data_set]['EclipticLongitude'][index]-std_list[data_set]['EclipticLongitude'][index],mean_list[data_set]['EclipticLongitude'][index]+std_list[data_set]['EclipticLongitude'][index])

# area_list = []
# for data_set in range(len(save_names)):
#     area = []
#     for i in range(len(mean_list[data_set]['EclipticLongitude'])):
#         if i % int(len(mean_list[data_set]['EclipticLongitude'])/10) == 0:
#             print(np.round(i/len(mean_list[data_set]['EclipticLongitude']),2))
#         area.append(integrate.nquad(f, [[np.arcsin(mean_list[data_set]['EclipticLatitude'][i]-std_list[data_set]['EclipticLatitude'][i]), np.arcsin(mean_list[data_set]['EclipticLatitude'][i]+std_list[data_set]['EclipticLatitude'][i])],
#                                          [mean_list[data_set]['EclipticLatitude'][i]-std_list[data_set]['EclipticLatitude'][i], mean_list[data_set]['EclipticLatitude'][i]+std_list[data_set]['EclipticLatitude'][i]]],
#                                         ))
#     area_list.append(area)
# pickle.dump(area_list, open(SAVEPATH+'/area_list' +save_name+end_string+'.npy', "wb"))


parameter = 'EclipticLongitude'
parameter2 = 'EclipticLatitude'
fig = plt.figure(figsize=fig_size)

for data_set in range(len(save_names)):
    plt.loglog(found_sources_list[data_set]['Amplitude'], std_list[data_set][parameter],  '.', label=labels_plot[data_set])
plt.xlabel('Amplitude')
plt.ylabel('Standard deviation ' +parameter)
plt.tight_layout()
plt.show()


for data_set in range(len(save_names)):
    for parameter in parameters:
        std_list[data_set][parameter] = np.array(std_list[data_set][parameter])
        mean_list[data_set][parameter] = np.array(mean_list[data_set][parameter])


for data_set in range(len(save_names)):
    std_list[data_set]['EclipticLatitude_no_sin'] = (np.arcsin(std_list[data_set]['EclipticLatitude']+ mean_list[data_set]['EclipticLatitude']) - np.arcsin(mean_list[data_set]['EclipticLatitude'])) * 2
# parameter = 'EclipticLongitude'


for data_set in range(len(save_names)):
    std_list[data_set]['angular_distance'] = angular_distance(mean_list[data_set]['EclipticLongitude']-std_list[data_set]['EclipticLongitude'], np.arcsin(mean_list[data_set]['EclipticLatitude']-std_list[data_set]['EclipticLatitude']), mean_list[data_set]['EclipticLongitude']+std_list[data_set]['EclipticLongitude'], np.arcsin(mean_list[data_set]['EclipticLatitude']+std_list[data_set]['EclipticLatitude']))


for data_set in range(len(save_names)):
    std_list[data_set]['angular_distance_latitude'] = angular_distance(mean_list[data_set]['EclipticLongitude'], np.arcsin(mean_list[data_set]['EclipticLatitude']-std_list[data_set]['EclipticLatitude']), mean_list[data_set]['EclipticLongitude'], np.arcsin(mean_list[data_set]['EclipticLatitude']+std_list[data_set]['EclipticLatitude']))


for data_set in range(len(save_names)):
    std_list[data_set]['angular_distance_longitude'] = angular_distance(mean_list[data_set]['EclipticLongitude']-std_list[data_set]['EclipticLongitude'], np.arcsin(mean_list[data_set]['EclipticLatitude']), mean_list[data_set]['EclipticLongitude']+std_list[data_set]['EclipticLongitude'], np.arcsin(mean_list[data_set]['EclipticLatitude']))

# parameter = 'EclipticLatitude_no_sin'
parameter = 'angular_distance_latitude'
n_bins = 30
fig = plt.figure(figsize=fig_size)

for data_set in range(len(save_names)):
    plt.hist(std_list[data_set][parameter] * 180/np.pi, bins=np.logspace(-2,2, n_bins),  density=False, histtype='step')
plt.xscale('log')
plt.xlabel(parameter)
plt.savefig(SAVEPATH+'/Evaluation/'+parameter+'std_histogram'+save_names[data_set]+end_string)
plt.tight_layout()
plt.show()

# area_list = pickle.load(open(SAVEPATH+'/area_list' +save_name+end_string+'.npy', "rb"))
# parameter = 'EclipticLongitude'
# parameter = 'EclipticLatitude'
# n_bins = 30
# fig = plt.figure(figsize=fig_size)
# for data_set in range(len(save_names)):
#     plt.hist(np.array(area_list[data_set])[:,0]*(180/np.pi)**2, bins=np.logspace(-3,3, n_bins), label=labels_plot[data_set], color=colors[data_set],  linestyle=linestyle[data_set], density=False, histtype='step', linewidth=line_width)
# plt.xscale('log')
# plt.xlim(10**-3,10**3)
# plt.xlabel('$\sigma_\mathrm{area} (\mathrm{deg}^2)$')
# plt.ylabel('Count')
# plt.grid(True)
# plt.legend(custom_lines, labels_plot, loc='upper left')
# plt.tight_layout()
# plt.savefig(SAVEPATH+'/Evaluation/confidence_area_std_histogram'+save_names[data_set]+end_string)
# plt.show()


# X = found_sources['Frequency']
# # mean['Frequency']
# Y = found_sources['Amplitude']
# Z = np.array(std[parameter])# * (180/np.pi)

# fig, ax = plt.subplots(1,1, figsize=fig_size)
# im = ax.scatter(X,Y,c=Z, norm= LogNorm(),  s=3)
# # fig.colorbar(im, ax=ax, label='std ' + labels[parameter])
# fig.colorbar(im, ax=ax, label='std ' + labels[parameter])
# # ax.set_xlim(0,2*np.pi)
# # ax.set_ylim(-np.pi/2,np.pi/2)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel('$f$')
# ax.set_ylabel(labels['Amplitude'])
# plt.savefig(SAVEPATH+'/Evaluation/'+parameter+'std'+save_names[data_set]+end_string)
# plt.show()

for data_set in range(len(save_names)):
    # pGB_injected_matched_list[data_set]['EclipticLatitude'] = np.sin(pGB_injected_matched_list[data_set]['EclipticLatitude'])
    for i in range(len(pGB_injected_matched_list[data_set]['EclipticLongitude'])):
        if pGB_injected_matched_list[data_set]['EclipticLongitude'][i] > np.pi:
            pGB_injected_matched_list[data_set]['EclipticLongitude'][i] -= 2*np.pi
        # pGB_injected_matched_list[data_set]['EclipticLongitude'][i] = pGB_injected_matched_list[data_set]['EclipticLongitude'][i]

found_sources_list = []
mean_list = []
std_list = []
in_contour_list = []

for data_set in range(len(save_names)):
    onlyfiles = [f for f in os.listdir(chain_paths[data_set]) if os.path.isfile(os.path.join(chain_paths[data_set], f))]
    found_sources = found_sources_matched_list[data_set]
    # found_sources = found_sources.sort_values(by=['Frequency'])
    # found_sources = found_sources.reset_index()
    # pGB_injected = pGB_injected_matched_list[data_set]
    # pGB_injected = pGB_injected.sort_values(by=['Frequency'])
    # pGB_injected = pGB_injected.reset_index()
    std = {}
    mean = {}
    in_contour = {}
    for parameter_x, parameter_y in itertools.combinations(parameters, 2):
        in_contour[parameter_x,parameter_y] = []
    for parameter in parameters:
        mean[parameter] = []
        std[parameter] = []
        in_contour[parameter] = []
        in_contour[parameter, parameter] = []
    for index in range(len(found_sources)):
        # if index > 46:
        #     continue
        # index_of_signal = np.searchsorted(found_sources['Frequency'], frequencies_in_folder[index])
        # if index_of_signal == len(onlyfiles): 
        #     index_of_signal = index_of_signal-1
        # if index_of_signal < len(onlyfiles)-1: 
        #     if abs(found_sources['Frequency'][index_of_signal]-frequencies_in_folder[index]) > abs(found_sources['Frequency'][index_of_signal+1]-frequencies_in_folder[index]):
        #         index_of_signal = index_of_signal+1
        # if index_of_signal != 0: 
        #     if abs(found_sources['Frequency'][index_of_signal]-frequencies_in_folder[index]) > abs(found_sources['Frequency'][index_of_signal-1]-frequencies_in_folder[index]):
        #         index_of_signal = index_of_signal-1
        index_of_signal = index
        # index_of_signal = 4804
        # index_of_signal = 502
        # index_of_signal = 9
        if index % int(len(onlyfiles)/10) == 0:
            print(np.round(index/len(onlyfiles),2))
        # print('frequency'+str(int(np.round(found_sources['Frequency'][index_of_signal]*10**9)))+'Sangria_1_full_cut.csv')
        try:
            df = pd.read_csv(chain_paths[data_set]+'frequency'+str(int(np.round(found_sources['Frequency'][index_of_signal]*10**12)))+'pHz'+ save_names[data_set] +'.csv')
        except:
            pass
        
        df['Frequency'] *= 10**-3
        # df['EclipticLatitude'] = np.sin(df['EclipticLatitude'])

        # for parameter_x, parameter_y in itertools.combinations(parameters, 2):
        #     in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region_2d_kde(parameter_x, parameter_y))
        # for parameter_x, parameter_y in itertools.combinations(parameters, 2):
        #     in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region_2D(parameter_x, parameter_y))
        # for parameter_x in parameters:
        #     in_contour[parameter_x,parameter_x].append(check_if_in_confidence_region_2d_kde(parameter_x, parameter_x))
        for parameter in parameters:
            in_contour[parameter].append(check_if_in_confidence_region_1D_std(parameter))
            
        # parameter_x = 'Frequency'
        # parameter_y = 'FrequencyDerivative'
        # in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region(parameter_x, parameter_y))
        # parameter_x = 'Frequency'
        # parameter_y = 'EclipticLongitude'
        # in_contour_f_ip_list.append(check_if_in_confidence_region(parameter_x, parameter_y))
        # parameter_y = 'EclipticLongitude'
        # parameter_x = 'EclipticLatitude'
        # in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region(parameter_x, parameter_y))
        # parameter_x = 'Amplitude'
        # parameter_y = 'Inclination'
        # in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region(parameter_x, parameter_y))

        # if in_contour_f_fd_list[-1] and in_contour_sky_list[-1] and in_contour_amp_inclination_list[-1]:
        #     in_contour_list.append(True)
        # else:
        #     in_contour_list.append(False)
    in_contour_list.append(in_contour)

# for index in range(len(onlyfiles)):
#     # if index > 46:
#     #     continue
#     index_of_signal = np.searchsorted(found_sources_matched_list[data_set]['Frequency'], frequencies_in_folder[index])
#     if abs(found_sources_matched_list[data_set]['Frequency'][index_of_signal]-frequencies_in_folder[index]) > abs(found_sources_matched_list[data_set]['Frequency'][index_of_signal+1]-frequencies_in_folder[index]):
#         index_of_signal = index_of_signal+1
#     elif abs(found_sources_matched_list[data_set]['Frequency'][index_of_signal]-frequencies_in_folder[index]) > abs(found_sources_matched_list[data_set]['Frequency'][index_of_signal-1]-frequencies_in_folder[index]):
#         index_of_signal = index_of_signal-1

#     # index_of_signal = 4804
#     # index_of_signal = 502
#     # index_of_signal = 9
#     if index % int(len(onlyfiles)/10) == 0:
#         print(np.round(index/len(onlyfiles),2))
#     # print('frequency'+str(int(np.round(found_sources_matched_list[data_set]['Frequency'][index_of_signal]*10**9)))+'Sangria_1_full_cut.csv')
#     try:
#         df = pd.read_csv(chain_path+'frequency'+str(int(np.round(found_sources_matched_list[data_set]['Frequency'][index_of_signal]*10**12)))+'pHzSangria_1year_dynamic_noise.csv')
#     except:
#         pass
#     injected_parameters = pGB_injected_matched_list[data_set]
#     if injected_parameters['EclipticLongitude'][index_of_signal] > np.pi:
#         injected_parameters.loc[index_of_signal, 'EclipticLongitude'] -= 2*np.pi

#     df['Frequency'] *= 10**-3


    # for parameter_x, parameter_y in itertools.combinations(parameters, 2):
    #     in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region(parameter_x, parameter_y))
    # for parameter_x, parameter_y in itertools.combinations(parameters, 2):
    #     in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region_2D(parameter_x, parameter_y))
    # for parameter_x in parameters:
    #     in_contour[parameter_x,parameter_x].append(check_if_in_confidence_region(parameter_x, parameter_x))
    # for parameter in parameters:
    #     in_contour[parameter].append(check_if_in_confidence_region_1D_kde(parameter))
        
    # parameter_x = 'Frequency'
    # parameter_y = 'FrequencyDerivative'
    # in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region(parameter_x, parameter_y))
    # parameter_x = 'Frequency'
    # parameter_y = 'EclipticLongitude'
    # in_contour_f_ip_list.append(check_if_in_confidence_region(parameter_x, parameter_y))
    # parameter_y = 'EclipticLongitude'
    # parameter_x = 'EclipticLatitude'
    # in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region(parameter_x, parameter_y))
    # parameter_x = 'Amplitude'
    # parameter_y = 'Inclination'
    # in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region(parameter_x, parameter_y))

    # if in_contour_f_fd_list[-1] and in_contour_sky_list[-1] and in_contour_amp_inclination_list[-1]:
    #     in_contour_list.append(True)
    # else:
    #     in_contour_list.append(False)


for data_set in range(len(save_names)):
    print('confidence region', confidence_threshold)
    for parameter_x, parameter_y in itertools.combinations(parameters, 2):
        print(parameter_x, parameter_y, int(np.round(len(np.ones_like(in_contour_list[data_set][parameter_x,parameter_y])[in_contour_list[data_set][parameter_x,parameter_y]])/len(in_contour_list[data_set][parameter_x,parameter_y])*100)), '%')

for data_set in range(len(save_names)):
    print(save_names[data_set])
    for parameter_x in parameters:
        # print(parameter_x, parameter_x, len(np.ones_like(in_contour[parameter_x,parameter_x])[in_contour[parameter_x,parameter_x]])/len(in_contour[parameter_x,parameter_x]))
        print(parameter_x, int(np.round(len(np.ones_like(in_contour_list[data_set][parameter_x])[in_contour_list[data_set][parameter_x]])/len(in_contour_list[data_set][parameter_x])*100)), '%')

# f_trues = len(np.ones_like(in_contour_f_fd_list)[in_contour_f_fd_list])
# sky_trues = len(np.ones_like(in_contour_sky_list)[in_contour_sky_list])
# amp_incl_trues = len(np.ones_like(in_contour_amp_inclination_list)[in_contour_amp_inclination_list])
# all_trues = len(np.ones_like(in_contour_list)[in_contour_list])

# print('frequency and derivative inside', f_trues/len(in_contour_f_fd_list))
# print('skylocation inside', sky_trues/len(in_contour_f_fd_list))
# print('amplitude and inclination inside', amp_incl_trues/len(in_contour_f_fd_list))
# print('all inside', all_trues/len(in_contour_f_fd_list))
# pdf = fastKDE.pdf_at_points(normalize(df[parameter_x]),normalize(df[parameter_y]),list_of_points = [(0.5,0.5)], axes=[ax,ax])



print('end')