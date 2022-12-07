# Tutorials

## Data quicklook

We provide basic tutorials which illustrate how to load the data for
the different LDC, and have a first quicklook. The data model differs
from one challenge to another to adapt to the specific needs of each
LDC, and also to follow the software tools evolution and new
conventions adopted within the consortium. 

For each of them, one will find instructions on how to load the data:

- using the standard python tools (h5py, numpy, ...)
- using the LDC toolbox.

- LDC1-Radler:
     - [standard python tools only](https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/notebooks/LDC1a-Radler.ipynb)
	 - [using LDC toolbox]() **TODO Maude**

- LDC2a-Sangria: 
     - [standard python tools only](https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/notebooks/LDC2a-Sangria.ipynb)
	 - [using LDC toolbox](https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/notebooks/sangria-1st-steps.ipynb)

- LDC2b-Spritz: 
     - [standard python tools only](https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/notebooks/LDC2b-Spritz.ipynb)

## Data analysis

### LDC1a Radler

Tutorials to get familiar with LISA data analysis have been produced
in the context of the first LISA Data Challenge (LDC1-Radler).

Those tutorials are based on the LDC1-Radler dataset and
[MLDC](https://gitlab.in2p3.fr/stas/MLDC) software, which was used to
produce the data.

Notebooks, hosted on Google colab, contain both data and software, and
can direclty be played by anyone:

- [search for verification binary (new)](https://colab.research.google.com/drive/1urEIq4nF9XO6qVZvv3BhowZoiCv2XFNA)
- [search for verification binary (old)](https://colab.research.google.com/drive/1-wor0tPPv_6GBYBB8BVZIirAKUR9k9hp)
- [search for MBHB](https://colab.research.google.com/drive/1jkBlgb_T3oWmozpiBmSARqs6pJ93ekpY)

### LDC2

On LDC2a-Sangria data:
- [search for one gb](https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/notebooks/search_for_one_gb.ipynb)
- [search for mbhb](https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/notebooks/mbhb_sangria_tutorial.ipynb) 

## Data generation

Data generation tools are developed along the different data
challenges, with increasing complexity and new features.

LDC1-Radler was produced with
[MLDC](https://gitlab.in2p3.fr/stas/MLDC) software for GW part and
[LISACode](https://gitlab.in2p3.fr/elisadpc/LISACode) for noise part.

LDC2a-Sangria was produced with the new
[LDC](https://gitlab.in2p3.fr/LISA/LDC) for GW part and
[LISANode](https://gitlab.in2p3.fr/j2b.bayle/LISANode) for noise part.

LDC2b-Spritz will be produced with the same LDC/LISANode scheme but
with a new interface which makes use of the latest tools and features
of the simulation toolkit:

- [lisa glitch](https://gitlab.in2p3.fr/lisa-simulation/glitch)
- [lisa orbits](https://gitlab.in2p3.fr/lisa-simulation/orbits)
- [lisa constants](https://gitlab.in2p3.fr/lisa-simulation/constants)

Future LDC data production will make use of additional tools, in
particular those developed within the INREP framework, like
[pyTDI](https://gitlab.in2p3.fr/LISA/LDPG/wg6_inrep/pytdi).

We maintain a
[notebook](https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/notebooks/gw-in-lisanode.ipynb)
with the latest functionalities which illustrates the different ways
of simulating the data (**TODO Maude update it**). As it makes use of
the latest features developed by the collaboration, this tutorial
requires a gitlab access to all software, given by the
[LISA full member](https://gitlab.in2p3.fr/LISA/full-members) gitlab
group membership.

To go further, refer to:
- [LDC tutorial session at LDPG 2021 workshop (slides+recording)](https://confluence.cnes.fr/display/LISALDPG/LDPG+Workshop+28-29+June%2C+1rst-2nd+July)
- [the data generation documentation](https://lisa.pages.in2p3.fr/LDC/data_generation)

