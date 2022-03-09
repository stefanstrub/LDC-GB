import numpy as np
import os
from os import listdir
from os.path import isfile, join

# get current directory
path = os.getcwd()
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)
folderpath = grandparent+"/LDC/pictures/LDC1-4s3"

onlyfiles = [f for f in listdir(folderpath) if isfile(join(folderpath, f))]
save_name = 'LDC1-4 low frequency'
search_range = [0.0003, 0.033333333333]
found_sources_mp_even_unsorted = []
frequencies = []
for i in range(len(onlyfiles)):
    print(i)
    sources = np.load(folderpath+'/'+onlyfiles[i], allow_pickle = True).tolist()
    frequency = 0
    j = 0
    while frequency == 0:
        try:
            frequency = sources[j][1][0][0][0]['Frequency']
        except:
            pass
        j += 1
    frequencies.append(frequency)
    found_sources_mp_even_unsorted.append(np.load(folderpath+'/'+onlyfiles[i], allow_pickle = True).tolist())
sorted_indexes = np.argsort(frequencies)
found_sources_mp_even = []
for i in range(len(onlyfiles)):
    found_sources_mp_even += found_sources_mp_even_unsorted[sorted_indexes[i]]

np.save(folderpath+'/found_sources_even3s'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', found_sources_mp_even)
