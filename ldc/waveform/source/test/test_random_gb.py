from ldc.waveform.source import SourceMaker
import numpy as np
import matplotlib.pyplot as plt
import sys



## Detached Galactic Binaries
## with gaussian randomization of mass, freq, incl
cfg_dgb = {'approximant': 'TD_fdot',
          'catalogs': '007_SeBa_r105_ag_wdwd_pop_highres_P025g70_N1e1.npy',
          'dfdt_GW' : True,
          'random_inclination': "gaussian_0.09",
          'random_mass' : "gaussian_0.01percent",
          'random_frequency' : "gaussian_0.01percent",
          'nsource': 2,
          'seed': 1234,
          'source_type': 'GB'}

## Interacting Galactic Binaries
## with gaussian randomization of mass and freq
cfg_igb = {'approximant': 'TD_fdot',
          'catalogs': 'AMCVn_GWR_MLDC_bulgefix_opt_N1e1.npy',
          'dfdt_GW' : False,
          'random_inclination': "uniform",
          'random_mass' : "gaussian_0.01percent",
          'random_frequency' : "gaussian_0.01percent",
          'nsource': 2,
          'seed': 1234,
          'source_type': 'GB'}


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--check-distribution', action='store_true',
                        help="Check and plot distribution of parameters")
    args = parser.parse_args()

    # Choose sources from catalogs with randomization
    for cfg in [cfg_dgb, cfg_igb]:
        source_maker = SourceMaker.type(cfg["source_type"],
                                        cfg["approximant"],
                                        catalogs=[cfg["catalogs"]]) 
        cat = source_maker.draw_random_catalog(**cfg)
        
        source_maker.close_logger()
        print(cat)
        print(cat.dtype)

    if args.check_distribution:

        for cfg, key in zip([cfg_dgb, cfg_igb], ["detached", "interacting"]):
        
            cfg['nsource'] = 1
            cfg_ref = cfg.copy() # for reference (no randomization)
            cfg_ref['random_mass'] = "0"
            cfg_ref['random_frequency'] = "0"
            if key == "detached":
                cfg_ref['random_inclination'] = "0"

            source_maker = SourceMaker.type(cfg_ref["source_type"],
                                            cfg_ref["approximant"],
                                            catalogs=[cfg_ref["catalogs"]]) 
            catref = source_maker.draw_random_catalog(**cfg_ref)
            source_maker.close_logger()

            # Check distribution DWD GB 
            Nrand = 5000
            cats = np.zeros(Nrand, dtype=cat.dtype)
            source_maker = SourceMaker.type(cfg["source_type"],
                                            cfg["approximant"],
                                            catalogs=[cfg["catalogs"]], verbose=False) 
            for i in range(Nrand):
                cfg['seed'] = cfg['seed'] + 1
                cat = source_maker.draw_random_catalog(**cfg)
                cats[i] = cat

            source_maker.close_logger()

            for xR in catref.dtype.names[1:]: 
                plt.figure()
                plt.hist(cats[xR], bins=50)
                plt.title("%sGB: "%(key) + xR + ": " + str(catref[xR][0]))

