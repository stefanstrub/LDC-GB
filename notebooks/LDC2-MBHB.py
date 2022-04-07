import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import copy
from tqdm import tqdm as tqdm
import xarray as xr
from astropy import units as u
import pandas as pd
import h5py

from ldc.common.constants import ldc_cosmo as cosmo
from ldc.lisa.orbits import Orbits
from ldc.lisa.projection import ProjectedStrain


plt.style.use(['seaborn-ticks','seaborn-deep'])
# %pylab inline
mpl.rcParams.update({'font.size': 16})
plt.rcParams['axes.grid'] = True
plt.rcParams["figure.figsize"] = (12,7)

from ldc.common.series import TimeSeries, FrequencySeries, TDI
from ldc.common.tools import compute_tdi_snr, window

sangria_fn = "LDC2_sangria_training_v2.h5"
tdi_ts = TDI.load(sangria_fn, name="obs/tdi")
dt = tdi_ts["X"].attrs["dt"]

# Build frequencyseries object for X,Y,Z
tdi_fs = TDI(dict([(k,tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

