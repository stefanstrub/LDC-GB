import numpy as np
import time
import multiprocessing as mp

# Prepare data
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[200000, 5])
data = arr.tolist()
print(data[:5])

# Solution Without Paralleization

def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
        time.sleep(0.0000000005) 
    return count

# results = []
# start = time.time()
# for row in data:
#     results.append(howmany_within_range(row, minimum=4, maximum=8))

# print(results[:10],time.time()-start)

# Parallelizing using Pool.apply()


# start = time.time()
# # Step 1: Init multiprocessing.Pool()
# pool = mp.Pool(mp.cpu_count()-2)

# # Step 2: `pool.apply` the `howmany_within_range()`
# results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]

# # Step 3: Don't forget to close
# pool.close()    

# print(results[:10],time.time()-start)

# Parallelizing using Pool.map()

# Redefine, with only 1 mandatory argument.
def howmany_within_range_rowonly(row, minimum=4, maximum=8):
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
        time.sleep(0.0000000005) 
    return count


start = time.time()
pool = mp.Pool(mp.cpu_count())

results = pool.map(howmany_within_range_rowonly, [row for row in data])

pool.close()

print(results[:10],time.time()-start)

# Parallelizing with Pool.starmap()
start = time.time()
pool = mp.Pool(mp.cpu_count())

results = pool.starmap(howmany_within_range, [(row, 4, 8) for row in data])

pool.close()

print(results[:10],time.time()-start)

# Parallel processing with Pool.apply_async()
start = time.time()
pool = mp.Pool(mp.cpu_count())

results = []

# Step 1: Redefine, to accept `i`, the iteration number
def howmany_within_range2(i, row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
        time.sleep(0.0000000005) 
    return (i, count)


# Step 2: Define callback function to collect the output in `results`
def collect_result(result):
    global results
    results.append(result)


# # Step 3: Use loop to parallelize
# for i, row in enumerate(data):
#     pool.apply_async(howmany_within_range2, args=(i, row, 4, 8), callback=collect_result)

# # Step 4: Close Pool and let all the processes complete    
# pool.close()
# pool.join()  # postpones the execution of next line of code until all processes in the queue are done.

# # Step 5: Sort results [OPTIONAL]
# results.sort(key=lambda x: x[0])
# results_final = [r for i, r in results]

# print(results_final[:10],time.time()-start)

# Parallel processing with Pool.apply_async() without callback function
# pool = mp.Pool(mp.cpu_count())

# results = []

# start = time.time()
# # call apply_async() without callback
# result_objects = [pool.apply_async(howmany_within_range2, args=(i, row, 4, 8)) for i, row in enumerate(data)]

# # result_objects is a list of pool.ApplyResult objects
# results = [r.get()[1] for r in result_objects]

# pool.close()
# pool.join()
# print(results[:10],time.time()-start)

# Parallelizing with Pool.starmap_async()

pool = mp.Pool(mp.cpu_count())

results = []
start = time.time()

results = pool.starmap_async(howmany_within_range2, [(i, row, 4, 8) for i, row in enumerate(data)]).get()

# With map, use `howmany_within_range_rowonly` instead
# results = pool.map_async(howmany_within_range_rowonly, [row for row in data]).get()

pool.close()
print(results[:10],time.time()-start)