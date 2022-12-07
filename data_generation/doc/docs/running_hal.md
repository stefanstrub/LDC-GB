# How to run a pipeline on the CNES computing center

The CNES documentation is available
[here](https://gitlab.cnes.fr/hpc/wikiHPC/wikis/).

## Setting the proxy 

Follow the
[instructions](https://gitlab.cnes.fr/hpc/wikiHPC/wikis/connexion-proxy#astuce)
to set a script for setting the `http_proxy` variables. Be careful to
replace the `_user` variable with your proxy username:
`_user='http-lisa-username'`.

Source this script if you need to access external resources (git
clone, wget, singularity pull): `source set_proxy.sh`

For git clone, do not forget to use the `https` url.

## Setting the environment

Please refer to the [Running at CCIN2P3](running_cc.md) documentation
to set a virtual environement to install snakemake.

Add the following lines to your `.bashrc` file: 
```
module load python/3.7.2
module load fftw
module load singularity
source venv/bin/activate
export SINGULARITY_BINDPATH="/work/SC/lisa/LDC/ancillary_data:/data"
```

## Running the pipeline

A template script `job_hal.sh` is provided.

You can use `screen` to keep the main process alive as long as the
jobs are running.



