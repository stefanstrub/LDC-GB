import numpy as np
import h5py
from ldc.waveform.source import load_gb_catalog
import ldc.io.hdf5 as h5io

old_fn = "../spritz/data/VGB.h5"
new_fn = "../spritz/data/VGB_v2.h5"

cat, units = load_gb_catalog(old_fn)

# fix non ascii character and amplitude
i = np.where(cat["Name"]=='CD–30o11223')[0][0]
cat[i]["Name"] = 'CD-30o11223'
assert cat[i]["Amplitude"] == 4.150747034e-21
cat[i]["Amplitude"] = 4.150747034e-22

# fix fdot
new_fdot = dict({"HM_Cnc":(7.485285541230417e-16, 3.627364057349999e-16),
                 "V407_Vul":(2.78570005589649e-17, 9.777592578125e-18),
                 "SDSS_J065133.34+284423.4":(2.6272309341370955e-17, 1.6736657391071642e-17),
                 "ES_Cet":(5.53680665282078e-17, -1.6481168858455217e-17)})

for k,v in new_fdot.items():
    i = np.where(cat["Name"]==k)[0][0]
    assert cat[i]["FrequencyDerivative"] == v[0]
    cat[i]["FrequencyDerivative"] = v[1]
    

# remove None in units
for k,v in units.items():
    if v is None:
        units[k] = '1'

h5io.save_array(new_fn, cat, **units)
new_cat, descr = h5io.load_array(new_fn)
print(new_cat)
print(descr)
