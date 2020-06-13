import numpy as np
import ldc.io.npz as npzio
import ldc.io.hdf5 as hdfio
import ldc.io.yml as ymlio

# a numpy array with some attributes
t_min = 0
t_max = 1e4
dt = 5
info = "doc"
N = int(t_max/dt)
arr = np.random.randn(N)

npz_fn = "test.npz"
npzio.save_array(npz_fn, arr, name='strain', t_min=t_min, t_max=t_max, dt=dt, info=info)
print (npzio.load_array(npz_fn, name='strain'))

hdf5_fn = "test.h5"
hdfio.save_array(hdf5_fn, arr, name='strain', mode='w', t_min=t_min, t_max=t_max, dt=dt, info=info)
print (hdfio.load_array(hdf5_fn, name='strain'))

# a simple dictionary (not nested)
cfg = {'dt': 5.0,
       'initial_position': 0.0,
       'initial_rotation': 0.0,
       'nominal_arm_length': 2500000000.0,
       'orbit_type': 'analytic'}
yml_fn = "config.yml"
ymlio.save_config(yml_fn, cfg)


hdfio.save_config(hdf5_fn, cfg)
print(hdfio.load_config(hdf5_fn))

# a record array
cat = np.rec.fromarrays([np.arange(5), np.arange(5)], names=["a", "b"])

npzio.save_array(npz_fn, cat, name='catalog')
print(npzio.load_array(npz_fn, name='catalog'))

hdfio.save_array(hdf5_fn, cat, name='catalog')
print(hdfio.load_array(hdf5_fn, name='catalog'))


