# Installation

## Getting source files

Cloning the repository: `git clone git@gitlab.in2p3.fr:maudelejeune/LDCPipeline.git`

## Dependencies

The LDCPipeline relies on:

- the [LDC library](https://gitlab.in2p3.fr/LISA/LDC) (`master` branch)
- the [LISANode software](https://gitlab.in2p3.fr/j2b.bayle/LISANode) (`ldc-pipeline` branch)
- [snakemake](https://snakemake.readthedocs.io/en/stable/) (python >=3.6, snakemake > 5.7.4)

Some test scripts perform comparisons to previous implementations:

- the [MLDC software](https://gitlab.in2p3.fr/stas/MLDC.git) (`development` branch)
- the [LISACode software](https://gitlab.in2p3.fr/elisadpc/LISACode) (`master` branch) 

## Content of the repository

The LDCPipeline toolbox provides: 

- workflow python scripts to run composite pipelines, based on the
  `snakemake` workflow management tool (in `pipeline` directory).
- jupyter notebook tutorials (in `notebook` directory).
- standalone python scripts which can be used to run a specific box
  (in `script` directory)

To install those tools: 

- `pip install -r requirements.txt`
- `python setup.py install` 

## Using the docker image

TBD
