from setuptools import setup, find_namespace_packages, Command
from distutils.core import Extension
from Cython.Distutils import build_ext
import os

class build_liborbits(Command):
    def run(self):
        os.system("make -C ldc/lisa/orbits/lib")
        #os.system("make install")

try:
    import cython
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

if USE_CYTHON:
    sources = ["ldc/lisa/orbits/lib/pyorbits.pyx",
               "ldc/lisa/orbits/lib/orbits.cc",
               "ldc/lisa/orbits/lib/common.cc"]
else:
    sources = ["ldc/lisa/orbits/lib/pyorbits.cpp"]

import numpy
orbits_ext = Extension("_orbits",
                       sources=sources,
                       language="c++",
                       include_dirs=[numpy.get_include()],
                       extra_compile_args=["-std=gnu++11"])


setup(
    name='ldc',
    version='1',
    description='LISA Data Challenge software',
    long_description='TBD',
    author='ldc-dev',
    author_email='ldc-dev@lisamission.org',
    cmdclass={'build_ext': build_ext, "build_liborbits": build_liborbits}, 
    packages=find_namespace_packages(include=['ldc.lisa.*','ldc.common.*']),
    zip_safe=False,
    install_requires=['numpy'],
    ext_modules=[orbits_ext],

)

