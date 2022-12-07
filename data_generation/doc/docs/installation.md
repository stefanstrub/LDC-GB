# Installation

## Getting source files

Cloning the repository: `git clone git@gitlab.in2p3.fr:lisa/LDC.git`

## Dependencies

The LDC pipeline relies on:

- the [LDC library](https://gitlab.in2p3.fr/LISA/LDC) (`develop` branch) and dependencies
- the [LISANode software](https://gitlab.in2p3.fr/j2b.bayle/LISANode) (v1.3)
- the [lisa simulation tools](https://gitlab.in2p3.fr/lisa-simulation)
- [snakemake](https://snakemake.readthedocs.io/en/stable/) (python >=3.6, snakemake > 5.7.4)

## Content of the repository

The toolbox provides: 

- workflow python scripts to run composite pipelines, based on the
  `snakemake` workflow management tool 
- jupyter notebook tutorials 
- standalone python scripts which can be used to run a specific box
  (in `script` directory)

## Using the docker image

The LDC CI provides a docker image which contains all the
dependencies, tagged with the prefix `prod`:
`docker://gitlab-registry.in2p3.fr/lisa/ldc:prod-master`


