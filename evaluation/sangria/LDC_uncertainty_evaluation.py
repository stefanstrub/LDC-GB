# from matplotlib.lines import _LineStyle
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.font_manager
import scipy
import numpy as np
import time
import pandas as pd
import os
import h5py
import itertools
from KDEpy import FFTKDE

import corner
import numpy.lib.recfunctions as recf

from ldc.lisa.noise import AnalyticNoise
from gb_evaluation import GBEval

from astropy import units as u
import astropy.coordinates as coord


# from fastkde import fastKDE

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




parameters = [
    "Amplitude",
    "EclipticLatitude",
    "EclipticLongitude",
    "Frequency",
    "FrequencyDerivative",
    "Inclination"
]



def normalize(array, min, max):
    return (array-min)/(max-min)

def check_if_in_confidence_region_1D_kde(parameter_x):
    min_x = np.min(df[parameter_x])
    max_x = np.max(df[parameter_x])
    # Notice how bw (standard deviation), kernel, weights and grid points are set
    x, pdf = FFTKDE(bw=0.1, kernel='gaussian').fit(np.array(normalize(df[parameter_x], min_x, max_x)), weights=None).evaluate(2**8)

    hist_indexes_sorted = np.unravel_index(np.argsort(-pdf, axis=None), pdf.shape)
    pdf_sorted = pdf[hist_indexes_sorted]

    plt.figure()
    plt.plot(x, pdf); plt.tight_layout()
    plt.show()

    plt.figure()
    plt.hist(df[parameter_x].values)
    plt.show()

    sum = 0
    for i in range(len(pdf_sorted)):
        sum += pdf_sorted[i]
        if sum > confidence_threshold:
            largest_index_in_contour = i
            break

    contour = np.zeros_like(pdf)
    for i in range(largest_index_in_contour):
        contour[hist_indexes_sorted[0][i]] = 1
    true_parameter = pGB_injected_matched_flat_df[parameter_x][index_of_signal]
    index_true = np.searchsorted(x, true_parameter)

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
    true_parameter = pGB_injected_matched_flat_df[parameter_x][index_of_signal]
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
    true_parameter_x = pGB_injected_matched_flat_df[parameter_x][index_of_signal]
    true_parameter_y = pGB_injected_matched_flat_df[parameter_y][index_of_signal]
    index_true = [np.searchsorted(hist_values_axes[1], true_parameter_x), np.searchsorted(hist_values_axes[2], true_parameter_y)]

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

def check_if_in_confidence_region(pdf, truth, parameter_x, parameter_y):
    min_x = np.min(pdf[parameter_x])
    max_x = np.max(pdf[parameter_x])
    min_y = np.min(pdf[parameter_y])
    max_y = np.max(pdf[parameter_y])

    grid_points = 2**6

    # grid_points = 2**5+1
    # ax = np.linspace(-0.15,1.15,grid_points)
    # myppdf,axes = fastKDE.ppdf(normalize(pdf[parameter_x], min_x, max_x),normalize(pdf[parameter_y], min_y, max_y), axes=[ax,ax])
    # x = ax

    data = np.array([normalize(pdf[parameter_x], min_x, max_x),normalize(pdf[parameter_y], min_y, max_y)]).T
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
    x_normalized = normalize(truth[parameter_x], min_x, max_x)
    y_normalized = normalize(truth[parameter_y], min_y, max_y)
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


in_contour = {}
# parameter_pairs = itertools.combinations(parameters, 2)
parameter_pairs = [['Amplitude', 'Inclination'],['Frequency', 'FrequencyDerivative'],['EclipticLongitude', 'EclipticLatitude']]

for parameter_x, parameter_y in parameter_pairs:
    in_contour[parameter_x,parameter_y] = []

for parameter_x in parameters:
    in_contour[parameter_x,parameter_x] = []
for parameter in parameters:
    in_contour[parameter] = []

confidence_threshold = 0.68

workdir = "/home/stefan/LDC/Sangria/evaluation"
# gb = GBEval('apc-l2it', workdir, submitted_noise=True)
# gb = GBEval('msfc-montana', workdir, submitted_noise=True)
gb= GBEval('eth', workdir, submitted_noise=True)

gb.load_from_workspace()

names = ['Amplitude', 'Inclination', 'EclipticLatitude', 'EclipticLongitude',
            'Frequency', 'FrequencyDerivative']
for i_inj in range(len(gb.inj_cat)):
    if i_inj % 10 != 0:
        continue
    
    pdf = gb.get_pdf(i_inj=i_inj, full_chain=True)
    if pdf is None:
        continue
    pdf = pdf[names]

    truth = gb.inj_cat[i_inj]

    if i_inj % int(len(gb.inj_cat)/10) == 0:
        print(np.round(i_inj/len(gb.inj_cat),2))

    # if pGB_injected_matched_flat_df['EclipticLongitude'][i_inj_of_signal] > np.pi:
    #     pGB_injected_matched_flat_df.loc[i_inj_of_signal, 'EclipticLongitude'] -= 2*np.pi


    for parameter_x, parameter_y in parameter_pairs:
        in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region(pdf, truth, parameter_x, parameter_y))
    # for parameter_x, parameter_y in itertools.combinations(parameters, 2):
    #     in_contour[parameter_x,parameter_y].append(check_if_in_confidence_region_2D(parameter_x, parameter_y))
    # for parameter_x in parameters:
    #     in_contour[parameter_x,parameter_x].append(check_if_in_confidence_region(parameter_x, parameter_x))
    # for parameter in parameters:
    #     in_contour[parameter].append(check_if_in_confidence_region_1D_kde(parameter))

print('number of signals', len(in_contour['Frequency','FrequencyDerivative']))
print( 'confidence area:', int(confidence_threshold*100) , '%')
for parameter_x, parameter_y in parameter_pairs:
    print(parameter_x, parameter_y, int(np.round(len(np.ones_like(in_contour[parameter_x,parameter_y])[in_contour[parameter_x,parameter_y]])/len(in_contour[parameter_x,parameter_y]),2)*100) , '%')

breakpoint()