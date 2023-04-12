import numpy as np
import pandas as pd
import h5py

# Ints don't support NaN, so store size in float
param_names = {"generator": str, "size": np.float64, "dt": np.float64,
               "t0": np.float64, "t_inj": np.float64, "inj_point": str,
               "beta": np.float64, "level": np.float64}

# Visit function must return None to continue, so make class to store results
# Based on https://stackoverflow.com/questions/31146036/how-do-i-traverse-a-hdf5-file-using-h5py
class GlitchList:
    def __init__(self):
        self.dsets = []
    def __call__(self, name, dset):
        if "injection_count" not in dset.attrs:
            return
        glitches = pd.DataFrame(index=range(dset.attrs["injection_count"]),
            columns=param_names.keys())
        glitches = glitches.astype(param_names)
        for idx in range(dset.attrs["injection_count"]):
            for param in param_names.keys():
                glitches.loc[idx, param] = dset.attrs[f"inj{idx}_{param}"]
        self.dsets += [glitches]

if __name__ == "__main__":
    input_file = "Spritz/vgb-glitch.h5"
    glitches = GlitchList()
    with h5py.File(input_file, "r") as vgb:
        vgb.visititems(glitches)
    # Concatenate injection points to single df
    all_glitches = pd.concat(glitches.dsets, ignore_index=True)
    all_glitches.to_csv("Spritz/vgb-glitch-params.txt", sep="\t", index=False)
