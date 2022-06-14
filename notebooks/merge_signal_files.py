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
folderpath = grandparent+"/LDC/pictures/LDC1-4uneven"
folderpath = grandparent+"/LDC/pictures/LDC1-4_4mHz_Euler"
folderpath = grandparent+"/LDC/pictures/LDC1-4_half_even"

onlyfiles = [f for f in listdir(folderpath) if isfile(join(folderpath, f))]
save_name = 'LDC1-4 low frequency'
save_name = 'LDC1-4_half_even3'
search_range = [0.0003, 0.033333333333]
found_sources_mp_even_unsorted = []
frequencies = []
for i in range(len(onlyfiles)):
    print(i)
    sources = np.load(folderpath+'/'+onlyfiles[i], allow_pickle = True).tolist()
    frequency = 0
    for j in range(len(sources)):
        try:
            frequency = sources[j][1][0][0][0]['Frequency']
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

# np.save(folderpath+'/found_sources_even3s'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', found_sources_mp_even)

np.save(folderpath+'/found_sources' +save_name+'.npy', found_sources_mp)
