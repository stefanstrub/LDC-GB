import numpy as np
from ldc.io.npz import NPZ
from ldc.io.hdf5 import HDF5
from ldc.io.yml import YML

# a numpy array with some attributes
t_min = 0
t_max = 1e4
dt = 5
info = "doc"
N = int(t_max/dt)
arr = np.random.randn(N)

npz = NPZ("test.npz")
npz.save_array(arr, name='strain', t_min=t_min, t_max=t_max, dt=dt, info=info)
print (npz.load_array(name='strain'))

hdf5 = HDF5("test.h5", remove_if_exist=True)
hdf5.save_array(arr, name='strain', t_min=t_min, t_max=t_max, dt=dt, info=info)
print (hdf5.load_array(name='strain'))

# a simple dictionary (not nested)
cfg = {'dt': 5.0,
       'initial_position': 0.0,
       'initial_rotation': 0.0,
       'nominal_arm_length': 2500000000.0,
       'orbit_type': 'analytic'}
yml = YML("config.yml")
yml.save_config(cfg)


hdf5 = HDF5("test.h5", remove_if_exist=True)
hdf5.save_config(cfg)
print(hdf5.load_config())

# a record array
cat = np.rec.fromarrays([np.arange(5), np.arange(5)], names=["a", "b"])

npz.save_array(cat, name='catalog')
print(npz.load_array(name='catalog'))

hdf5.save_array(cat, name='catalog')
print(hdf5.load_array(name='catalog'))


