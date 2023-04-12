import numpy as np
import ldc.io.hdf5 as h5io

DATAPATH = "/home/maude/data/LDC/sangria/1.7"

for src in ["mbhb", "sum", "vgb", "dgb", "igb"]:
    data = []
    for name in ["X", "Y", "Z"]:
        X, attrs = h5io.load_array(DATAPATH+f"/{src}-tdi.h5", name=name) # or Y, Z
        ineg = np.where(X[:,0]>=0)[0][0] # remove negative time indices
        tvec = X[ineg:-1, 0]
        data.append(X[ineg:-1, 1])
    rec = np.rec.fromarrays([tvec]+data, names=["t", "X", "Y", "Z"])
    h5io.save_array(DATAPATH+f"/{src}-tdi-XYZ.h5", rec, name='tdi')
    
