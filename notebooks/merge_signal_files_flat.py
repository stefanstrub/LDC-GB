import numpy as np
import os
from os import listdir
from os.path import isfile, join
import pickle
import matplotlib.pyplot as plt
import sys

# get current directory
path = os.getcwd()
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4/found_signals"
folderpath_parent = grandparent+"/LDC/pictures/Sangria/found_signals"
# folderpath_parent = grandparent+"/LDC/Spritz/found_signals"
# folderpath_save = grandparent+"/LDC/pictures/LDC1-4"
folderpath_save = grandparent+"/LDC/pictures/Sangria"
# folderpath_save = grandparent+"/LDC/Spritz/"

which_run = str(sys.argv[1])+'_'
which_run = ''
week = int(sys.argv[5])
name = '_original_Sangria_'+str(week)+'w_mbhb_SNR9_'+which_run+'seed1'
save_name = 'Spritz' + name
# save_name = 'Sangria' + name
save_name =  name
folderpath = folderpath_parent + name
onlyfiles = [f for f in listdir(folderpath) if isfile(join(folderpath, f))]
found_sources_mp_even_unsorted = []
frequencies = []
times = []
for i in range(len(onlyfiles)):
    try:
        sources = pickle.load(open(folderpath+'/'+onlyfiles[i], 'rb'))
    except:
        sources = np.load(folderpath+'/'+onlyfiles[i], allow_pickle = True).tolist()
    frequency = 0
    for j in range(len(sources)):
        try:
            # frequency = sources[j][1][0][0][0]['Frequency']
            times.append(sources[j][5])
            frequency = sources[j][4][0]
            frequencies.append(frequency)
            found_sources_mp_even_unsorted.append(sources[j][3])

        except:
            # found_sources_mp_even_unsorted.append(sources[j][3])
            found_sources_mp_even_unsorted.append(sources[j])
    # frequencies.append(frequency)
    # found_sources_mp_even_unsorted.append(sources)
try:
    found_sources_in_flat = np.concatenate(found_sources_mp_even_unsorted)
except:
    found_sources_in_flat = found_sources_mp_even_unsorted
found_sources_in_flat_frequency = []
for i in range(len(found_sources_in_flat)):
    found_sources_in_flat_frequency.append(found_sources_in_flat[i]['Frequency'])
found_sources_in_flat_frequency = np.asarray(found_sources_in_flat_frequency)
found_sources_in_flat = np.asarray(found_sources_in_flat)
indexes_in = np.argsort(found_sources_in_flat_frequency)
found_sources_in_flat_frequency = found_sources_in_flat_frequency[indexes_in]
found_sources_in_flat = found_sources_in_flat[indexes_in]
print(len(found_sources_in_flat))

print('time', np.round(np.sum(times)/3600,0), 'hours ')
# print('max time', np.round(np.max(times)/60,1), 'minutes')

# plt.figure()
# plt.plot(found_sources_in_flat_frequency,found_sources_in_flat_frequency, '.')
# plt.show()


pickle.dump(found_sources_in_flat, open(folderpath_save+'/found_sources' +save_name+'_flat.pkl', 'wb'))
# np.save(folderpath_save+'/found_sources_' +save_name+'.npy', found_sources_mp)
