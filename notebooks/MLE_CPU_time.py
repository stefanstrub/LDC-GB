import os
import numpy as np


# sacct -a -j 15401519 --format=user%10,jobname%10,node%10,start%10,end%10,elapsed%10,CPUTime,MaxRSS
#144-194,196,197,198,199
#110, 114, 118, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 133, 134, 137, 138, 139, 140, 141, 143, 196, 197, 198, 199, 201, 202, 209, 213, 224
#224, 213, 209

#142,148,150,158,160,164,169

# get current directory
path = os.getcwd()
 
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)
Radler = False
if Radler:
    DATAPATH = grandparent+"/LDC/Radler/data/"
    SAVEPATH = grandparent+"/LDC/pictures/LDC1-4/duration/"
else:
    DATAPATH = grandparent+"/LDC/Sangria/data/"
    SAVEPATH = grandparent+"/LDC/pictures/Sangria/duration/"

# Define a function to convert time string to seconds
def time_to_seconds(time_str):
    # Split the time string into hours, minutes, and seconds
    hours, minutes, seconds = map(int, time_str.split(':'))
    # Calculate the total seconds
    total_seconds = hours * 3600 + minutes * 60 + seconds
    return total_seconds

time = 0
times = []
numer_of_times = 0
# Open the text file in read mode
with open(SAVEPATH+'times_sangria_12m_even.txt', "r") as file:
# with open(SAVEPATH+'out.txt', "r") as file:
    # Loop through each line in the file
    for line in file:
        # Split the line into a list of elements using whitespace as delimiter
        elements = line.split()
        # Check if the line has at least 10 elements
        if len(elements) >= 5:
            # Access the 10th element (index 9) and print it
            # print(time_to_seconds(elements[5]))
            numer_of_times += 1
            times.append(time_to_seconds(elements[5]))
            time += time_to_seconds(elements[5])
    time = np.sum(times)
    max_time = np.max(times)
    print('time', np.round(time/3600,1), 'hours')
    print('max time', np.round(max_time/60,1), 'minutes')
time /= 3
print(np.round(time/3600), 'CPU hours')
print(file.read())