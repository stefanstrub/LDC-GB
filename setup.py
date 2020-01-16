from setuptools import setup, Command
from distutils.core import Extension
from Cython.Distutils import build_ext
import os

try:
    from setuptools import find_namespace_packages
except ImportError:
    raise("setuptools>=41.0.0 is required")


class build_liborbits(Command):
    description = 'compile c++ orbits library'
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def sub_commands(self):
        pass
    def run(self):
        os.system("make -C ldc/lisa/orbits/lib")
        #os.system("make install")

try:
    import cython
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

if 1:#USE_CYTHON:
    orbits_sources = ["ldc/lisa/orbits/lib/pyorbits.pyx",
                      "ldc/lisa/orbits/lib/orbits.cc",
                      "ldc/lisa/orbits/lib/common.cc"]
else:
    orbits_sources = ["ldc/lisa/orbits/lib/pyorbits.cpp"]

import numpy
orbits_ext = Extension("_orbits",
                       sources=orbits_sources,
                       language="c++",
                       include_dirs=[numpy.get_include(), 'ldc/common/constants'],
                       extra_compile_args=["-std=gnu++11"])

fastGB_sources = ["ldc/waveform/fastGB/pyGB.pyx",
                  "ldc/waveform/fastGB/GB.c", "ldc/waveform/fastGB/LISA.c"]

fastGB_ext = Extension("fastGB",
                       sources=fastGB_sources,
                       language="c",
                       include_dirs=[numpy.get_include(), "ldc/common/constants"],
                       libraries=['fftw3', 'gsl'])#,

setup(
    name='ldc',
    version='1',
    description='LISA Data Challenge software',
    long_description='TBD',
    author='ldc-dev',
    author_email='ldc-dev@lisamission.org',
    cmdclass={'build_ext': build_ext, "build_liborbits": build_liborbits}, 
    packages=find_namespace_packages(include=['ldc.lisa.*','ldc.common.*', 'ldc.waveform.*']),
    zip_safe=False,
    install_requires=['numpy'],
    ext_modules=[orbits_ext, fastGB_ext],

)

