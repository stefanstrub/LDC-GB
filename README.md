# LISA Data Challenge software

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.7332221.svg)](https://doi.org/10.5281/zenodo.7332221)

## Installation of the latest released version

`pip install lisa-data-challenge`

## Installation of the dev version

### Cloning the gitlab project

The default working branch is named `develop`. 
`git clone -b develop https://gitlab.in2p3.fr/LISA/LDC.git`

### Prerequisites

Install the prerequisites :

- GSL : `apt-get install libgsl-dev` or `conda install gsl`
- FFTW3 : `apt-get install libfftw3-dev` or `conda install fftw`
- lisaconstants: `pip install git+https://gitlab.in2p3.fr/lisa-simulation/constants.git@master --no-binary :all:`
- few: see https://bhptoolkit.org/FastEMRIWaveforms/html/index.html

### Installation

Make sure that all requirements are met.

The `requirements.txt` file defines the reference version for most of
the dependencies for a python3.9 installation as recommended by
[LISA-CDE](https://gitlab.in2p3.fr/LISA/lisa-cde), but other versions
of the listed package might work. 

To comply with the CDE environement:
`pip install -r requirements.txt`

Paths to FFTW and GSL can be set explicitly by editing `setup.cfg`.

Extensions for specific fast waveform generator must be explicitely
listed in the installation command line. To have them all:

`python setup.py install --with-fastGB --with-imrphenomD --with-fastAK`


## Documentation

- [API documentation](https://lisa.pages.in2p3.fr/LDC/)

## Use policy

Do not forget to associate the authors of this software to your
research:

- Please cite the DOI (see badge above) and acknowledge the LDC
  working group in any publication which makes use of it

- Do not hesitate to send an email for help and/or collaboration:
  ldc-at-lisamission.org, ldc-chairs-at-lisamission.org
