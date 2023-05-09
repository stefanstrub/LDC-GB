import numpy as np
import os
from os import listdir
from os.path import isfile, join
import pickle

# get current directory
path = os.getcwd()
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)
folderpath_parent = grandparent+"/LDC/pictures/LDC1-4/found_signals"
# folderpath_parent = grandparent+"/LDC/pictures/Sangria/found_signals"
folderpath_save = grandparent+"/LDC/pictures/LDC1-4"
# folderpath_save = grandparent+"/LDC/pictures/Sangria"

name = '_24m'
save_name = 'Radler' + name
# save_name = 'Sangria' + name
folderpath = folderpath_parent + name
onlyfiles = [f for f in listdir(folderpath) if isfile(join(folderpath, f))]
found_sources_mp_even_unsorted = []
frequencies = []
for i in range(len(onlyfiles)):
    try:
        sources = pickle.load(open(folderpath+'/'+onlyfiles[i], 'rb'))
    except:
        sources = np.load(folderpath+'/'+onlyfiles[i], allow_pickle = True).tolist()
    frequency = 0
    for j in range(len(sources)):
        try:
            # frequency = sources[j][1][0][0][0]['Frequency']
            frequency = sources[j][4][0]
            frequencies.append(frequency)
            found_sources_mp_even_unsorted.append(sources[j])

        except:
            pass
    # frequencies.append(frequency)
    # found_sources_mp_even_unsorted.append(sources)
sorted_indexes = np.argsort(frequencies)
found_sources_mp = []
for i in range(len(sorted_indexes)):
    found_sources_mp.append(found_sources_mp_even_unsorted[sorted_indexes[i]])

try:
    time = 0
    times = []
    for i in range(len(found_sources_mp)):
        times.append(found_sources_mp[i][5])
    time = np.sum(times)
    max_time = np.max(times)
    print('time', np.round(time/3600,1), 'hours')
    print('max time', np.round(max_time/60,1), 'minutes')
except:
    print('no time')

np.save(folderpath_save+'/found_sources_' +save_name+'.npy', found_sources_mp)
