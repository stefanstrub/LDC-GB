from setuptools import setup, find_namespace_packages

setup(
    name='ldc',
    version='1',
    description='LISA Data Challenge software',
    long_description='TBD',
    author='ldc-dev',
    author_email='ldc-dev@lisamission.org',
    
    packages=find_namespace_packages(include=['ldc.lisa.*','ldc.common.*']),
    zip_safe=False,
)

