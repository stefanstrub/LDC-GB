from setuptools import setup, Command, find_namespace_packages
#from distutils.core import Extension
from setuptools import Extension
from Cython.Distutils import build_ext
import os
import subprocess
import distutils
import sys

GSL_CFLAGS = subprocess.check_output("gsl-config --cflags", shell=True)

try:
    import cython
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False
ext_c = '.pyx' if USE_CYTHON else '.c'
ext_cc = '.pyx' if USE_CYTHON else '.cpp'

import numpy
orbits_ext = Extension("_orbits",
                       sources=["ldc/lisa/orbits/lib/pyorbits.pyx",
                                "ldc/lisa/orbits/lib/orbits.cc",
                                "ldc/lisa/orbits/lib/common.cc"],
                       language="c++",
                       include_dirs=[numpy.get_include(), 'ldc/common/constants',
                                     'ldc/lisa/orbits/lib'],
                       extra_compile_args=["-std=gnu++11"])

fastGB_ext = Extension("fastGB",
                       sources=["ldc/waveform/fastGB/pyGB.pyx",
                                "ldc/waveform/fastGB/GB.cc",
                                "ldc/waveform/fastGB/LISA.cc",
                                "ldc/lisa/orbits/lib/c_wrapper.cc",
                                "ldc/lisa/orbits/lib/orbits.cc"], 
                       language="c++",
                       include_dirs=[numpy.get_include(), "ldc/common/constants",
                                     GSL_CFLAGS.decode()[2:],
                                     'ldc/lisa/orbits/lib'],
                       library_dirs=['ldc/lisa/orbits/lib'],
                       extra_compile_args = ['-fpermissive', "-std=gnu++11"],
                       #["-std=c99"],
                       libraries=['fftw3', 'gsl', 'gslcblas', 'orbits'])


imr_phenomd_ext = Extension("pyimrphenomD",
                            sources=["ldc/waveform/imrphenomD/pyimrphenomD"+ext_c,
                                     "ldc/waveform/imrphenomD/IMRPhenomD.c",
                                     "ldc/waveform/imrphenomD/IMRPhenomD_internals.c"],
                            include_dirs=[numpy.get_include(),
                                          'ldc/common/constants'],
                            language="c",
                            extra_compile_args=["-std=c99", "-O3"],
                            libraries=["gsl", "gslcblas"],
                            #library_dirs=[lib_gsl_dir])
                            )

fast_ak_ext = Extension("fastAK",
                        sources=["ldc/waveform/fastAK/fastAK"+ext_cc,
                                 "ldc/waveform/fastAK/EMRItemplate.cc",
                                 "ldc/waveform/fastAK/FreqAK_RAv2.cc",
                                 "ldc/lisa/orbits/lib/orbits.cc"],
                        include_dirs=[numpy.get_include(),
                                      'ldc/common/constants',
                                      'ldc/lisa/orbits/lib'],
                        library_dirs=['ldc/lisa/orbits/lib'],
                        language="c++",
                        extra_compile_args=["-std=gnu++11"],
                        libraries=["gsl", "gslcblas", 'orbits'])


default = [fastGB_ext, imr_phenomd_ext, fast_ak_ext]
ext_modules = [orbits_ext]
ext_modules+= [fastGB_ext, imr_phenomd_ext, fast_ak_ext]
for iopt, opt in enumerate(['--no-fastGB', '--no-imrphenomD', '--no-fastAK']):
    if opt in sys.argv:
        ext_modules.remove(default[iopt])
        sys.argv.remove(opt)

if orbits_ext in ext_modules:
    os.system("make -s -C ldc/lisa/orbits/lib")
    ext_modules.append(orbits_ext)

setup(cmdclass = {'build_ext': build_ext },
      ext_modules=ext_modules,
      packages=find_namespace_packages(include=['ldc.lisa.*','ldc.common.*',
                                                'ldc.waveform.*', 'ldc.io.*',
                                                'ldc.utils.*']),)
