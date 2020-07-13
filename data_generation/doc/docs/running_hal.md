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

Here is an example of a snakemake command line to run the pipeline on
HAL, using singularity container: 

`snakemake --use-singularity --cores 250 --default-resources "mem_mb=10000" --cluster "qsub -l select=1:ncpus=1:mem={resources.mem_mb}mb,walltime=24:00:00" --batch batch_merge=2/3`

You can use `screen` to keep the main process alive as long as the
jobs are running.

## Looking at the results with a jupyter notebook

At this stage, we suppose that you have install the ldc package into
your virtual environement.

Add a symbolic link from the lisa space to your home: `ln -s
/work/LC/lisa lisa`.

Set a dedicated kernel: `ipython kernel install --user --name ldc_kernel`

Tune the kernel configuration file
`.local/share/jupyter/kernels/ldc_kernel/kernel.json` to use a starter
script:
```
{
  "display_name": "ldc_kernel",
  "language":     "python",
  "argv": [
      "/home/u/username/jupyter-helper.sh",
      "-f",
      "{connection_file}"
  ]
}
```

Set this starter script `jupyter-helper.sh`: 
```
#!/bin/bash
module load python
module load fftw
source /home/u/username/venv/bin/activate
exec /home/u/username/venv/bin/python -m ipykernel_launcher "$@"
```

Go to https://jupyterhub.cnes.fr/ to connect and choose the ldc
kernel.



