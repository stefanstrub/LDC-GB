# LISA Data Challenge software

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.7332221.svg)](https://doi.org/10.5281/zenodo.7332221)

LDC provides a set of tools to generate and analyse the LDC datasets. 

## Installation of the latest released version

`pip install lisa-data-challenge`

## Installation of the dev version

### Cloning the gitlab project

The default working branch is named `develop`. 
`git clone -b develop https://gitlab.in2p3.fr/LISA/LDC.git`

### Installation

By default, `pyproject.toml` will be used to generate a temporary
environement to build the package (see
https://pip.pypa.io/en/stable/reference/build-system/pyproject-toml/)

`pip install .`

For older version of `setuptools`, a `setup.py` file is also provided.

`python setup.py install`

## Troubleshooting

### Prerequisites

- GSL : `apt-get install libgsl-dev` or `conda install gsl`
- FFTW3 : `apt-get install libfftw3-dev` or `conda install fftw`

Paths to FFTW and GSL can be set explicitly by editing `setup.cfg`.

### Python dependencies

Make sure that all requirements are met.

The `requirements.txt` file defines the reference version for most of
the dependencies for a python3.9 installation as recommended by
[LISA-CDE](https://gitlab.in2p3.fr/LISA/lisa-cde), but other versions
of the listed package might work. 

To comply with the CDE environement:
`pip install -r requirements.txt`

Extensions for specific fast waveform generator can be disabled in the
installation command line:

`python setup.py install --no-fastGB --no-imrphenomD --no-fastAK`

### Extra dependencies

Some external tools are interfaced by the LDC and need separate installation:

- EMRI waveform with few: see https://bhptoolkit.org/FastEMRIWaveforms/html/index.html
- Fast BH waveform with lisabeta: see https://gitlab.in2p3.fr/marsat/lisabeta

## Documentation

- [LISA Data Challenge portal](https://lisa-ldc.lal.in2p3.fr)
- [API documentation](https://lisa.pages.in2p3.fr/LDC/)

## Use policy

Do not forget to associate the authors of this software to your
research:

- Please cite the DOI (see badge above) and acknowledge the LDC
  working group in any publication which makes use of it

- Do not hesitate to send an email for help and/or collaboration:
  ldc-at-lisamission.org, ldc-chairs-at-lisamission.org
