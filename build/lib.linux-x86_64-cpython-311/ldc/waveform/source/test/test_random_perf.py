from ldc.waveform.source import SourceMaker
import numpy as np
import matplotlib.pyplot as plt
import sys
import time

## Verification Galactic Binairis
cfg_gb = {'approximant': 'TD_fdot',
          'catalogs': 'Catalog_NoID_single.npy',
          'nsource': 20000000,
          'seed': 1234,
          'source_type': 'GB'}

## Detached Galactic Binaries
## with gaussian randomization of mass, freq, incl
cfg_dgb = {'approximant': 'TD_fdot',
           'catalogs': '007_SeBa_r105_ag_wdwd_pop_highres_P025g70_8col.npy', 
           'dfdt_GW' : True,
           'random_inclination': "gaussian_0.09",
           'random_mass' : "gaussian_0.01percent",
           'random_frequency' : "gaussian_0.01percent",
           'nsource': 20000000,
           'seed': 1234,
           'source_type': 'GB'}

### Choose sources from catalogs without randomization
source_maker = SourceMaker.type(cfg_gb["source_type"],
                                cfg_gb["approximant"],
                                catalogs=[cfg_gb["catalogs"]]) 
t0 = time.time()
cat = source_maker.choose_from_catalog(**cfg_gb)
t1 = time.time()
source_maker.close_logger()

print("choose from cat (%d sources): %f sec"%(cfg_gb["nsource"], t1-t0))

source_maker = SourceMaker.type(cfg_dgb["source_type"],
                                cfg_dgb["approximant"],
                                catalogs=[cfg_dgb["catalogs"]], verbose=False) 
t0 = time.time()
cat = source_maker.draw_random_catalog(**cfg_dgb)
t1 = time.time()
source_maker.close_logger()

print("draw random cat (%d sources): %f sec"%(cfg_gb["nsource"], t1-t0))
 
