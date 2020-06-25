from setuptools import setup, Command
from distutils.core import Extension
from Cython.Distutils import build_ext
import os
import subprocess

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
                       include_dirs=[numpy.get_include(), 'ldc/common/constants'],
                       extra_compile_args=["-std=gnu++11"])

fastGB_ext = Extension("fastGB",
                       sources=["ldc/waveform/fastGB/pyGB.pyx",
                                "ldc/waveform/fastGB/GB.c",
                                "ldc/waveform/fastGB/LISA.c"], 
                       language="c",
                       include_dirs=[numpy.get_include(), "ldc/common/constants",
                                     GSL_CFLAGS.decode()[2:]],
                       extra_compile_args = ["-std=c99"], 
                       libraries=['fftw3', 'gsl'])#,

imr_phenomd_ext = Extension("pyimrphenomD",
                            sources=["ldc/waveform/imrphenomD/pyimrphenomD"+ext_c,
                                     "ldc/waveform/imrphenomD/IMRPhenomD.c",
                                     "ldc/waveform/imrphenomD/IMRPhenomD_internals.c"],
                            include_dirs = [numpy.get_include(), 'ldc/common/constants'],#, include_gsl_dir],
                            language="c",
                            extra_compile_args = ["-std=c99", "-O3"],
                            libraries=["gsl", "gslcblas"],
                            #library_dirs=[lib_gsl_dir])
                            )

setup(
    name='ldc',
    version='0.1',
    description='LISA Data Challenge software',
    long_description='TBD',
    author='ldc-dev',
    author_email='ldc-dev@lisamission.org',
    cmdclass={'build_ext': build_ext, "build_liborbits": build_liborbits}, 
    packages=find_namespace_packages(include=['ldc.lisa.*','ldc.common.*', 'ldc.waveform.*', 'ldc.io.*']),
    zip_safe=False,
    install_requires=['numpy'],
    ext_modules=[orbits_ext, fastGB_ext, imr_phenomd_ext],
    scripts=[os.path.join("data_generation/scripts", s) for s in ['arm_projection',
                                                                  "strain_combination",
                                                                  "strain_interpolation",
                                                                  'run_lisanode',
                                                                  'prep_lisanode',
                                                                  'source_selection',
                                                                  'data_release']],
    
)

