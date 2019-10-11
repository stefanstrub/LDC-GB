from setuptools import setup, find_namespace_packages
from distutils.core import Extension
from Cython.Distutils import build_ext

# TODO: if orbits + if cython, use precompiled .cpp otherwise
import numpy
orbits_ext = Extension("_orbits",
                       sources=["ldc/lisa/orbits/lib/pyorbits.pyx",
                                "ldc/lisa/orbits/lib/orbits.cc",
                                "ldc/lisa/orbits/lib/common.cc"],
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
    cmdclass={'build_ext': build_ext}, 
    packages=find_namespace_packages(include=['ldc.lisa.*','ldc.common.*']),
    zip_safe=False,
    install_requires=['numpy'],
    ext_modules=[orbits_ext],

)

