# Prerequisites

- GSL (e.g., `apt-get install libgsl-dev` or `conda install gsl`)
- FFTW3 (e.g., `apt-get install libfftw3-dev` or `conda install fftw`)

# Installation

- `pip install -r requirements.txt`
- `python setup.py build_liborbits install`
- `python setup.py install`

Paths to FFTW and GSL can be set explicitly by editing `setup.cfg`.

# Documentation

- [API documentation](https://lisa.pages.in2p3.fr/LDC/)

# Running the notebooks

To run the MLE search run the script called LDC1-3_MLE.py

To calculate the posterior distribution using Gaussian Process Regression run the script called LDC1-3_GP.py
