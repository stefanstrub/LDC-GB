# How to run a pipeline on the CCIN2P3 computing center

The CCIN2P3 documentation is available
[here](https://doc.cc.in2p3.fr).

## Setting the environment

In order to use snakemake, one has to install it, with a recent python
version. In the following, we will use `virtualenv` to set up a
working environment.

The default python3 version is 3.6.8

### Setting the virtualenv

```
> /usr/bin/pip3 install virtualenv --user
> cd /sps/hep/lisaf/username
> ~/.local/bin/virtualenv -p /usr/bin/python3 venv
> source /sps/hep/lisaf/username/venv/bin/activate
```

### Installing the pipeline dependancies

```
> pip install snakemake==5.7.4
```

### Setting up a new pipeline

The data are already available in the `/sps/hep/lisaf/mlejeune/data`
directory.

```
> cd pipeline
> runsim --init /sps/hep/lisaf/username/sangria
> cd /sps/hep/lisaf/username/sangria
> runsim --data /data
```

To download the singularity image: 

```
> ccenv singularity 3.5.2
> singularity pull --docker-username gitlab-deploy-token --docker-password ... ldcpipe.sif oras://gitlab-registry.in2p3.fr/lisa/ldc/ldc:prod-master
```

### Configuring the bash profile

In order to automatically activate this environment, put the following
lines in the `~/.bash_profile` file. 

```
source /sps/hep/lisaf/username/venv/bin/activate
ccenv singularity 3.5.2
export SINGULARITY_BINDPATH="/sps/hep/lisaf/mlejeune/data:/data"
```

### Using LISA project as default

Running the pipeline requires to use the `demon` queue, which is
restricted to some projects. It is thus important to use the LISA
accounting.

To know your default group, look at the `current group name` of
`newgroup --query`

To change group: `newgroup lisaf`


## Running the pipeline

### Galactic binaries

It is wise to run the galactic binaries subworkflow first, to adjust
the requested resources to its specifities (in particular, it requires
a very limited amount of memory and a high number of jobs, which is a
good configuration for the CCIN2P3 long queue).

The first step of the pipeline, which consists in choosing the source
in the 30 millions sources catalog should be run on the head node. 

```
cd gb-pipeline
snakemake -j 1 --use-singularity run1/dgb.npy
snakemake -j 1 --use-singularity run1/igb.npy
```

Then, a template script `job.sh` is provided. It can be tuned with
respect to the resources needed by the pipeline.

The first 9 lines should not be changed, they correspond to the
resources used by the main snakemake process. 

In the last line, set the `-j` option to the number of batch jobs to
want to run in parallel and the `h_rt` option to the expected walltime
of each job.

It is crutial to have the possibility to connect to a head node
without having to answer interactive question about host authenticity. 
Thus before running the job, run an `ssh cca004`. 

Once everything in place: `qsub job.sh`

#### Example

In this example, we include 1e6 sources. We limit the walltime of each
job to 6 hours. In 6 hours, around 4000 sources can be computed (5s
per source). Thus we will need to use around 250 batch jobs.

In `param_gb.yml` set `nsource` to 1e6 and in `config.yml` set
`nbatch` to 250. 

In order to run the pipeline in a single raw, set the `-j` option to
250, and the `h_rt` to '10:00:00'. 

You can also run this pipeline in several times, using the batch
option: 
`snakemake -j 50 --use-singularity --cluster "ssh username@cca004 /opt/sge/bin/lx-amd64/qsub -l h_rt=06:00:00,s_rss=1G -P P_lisaf" --batch batch_merge=1/5`

then 2/5, 3/5, 4/5 and 5/5.



