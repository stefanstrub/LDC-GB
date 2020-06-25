# Note on the performances and vectorization

This note explains the design choices which have been made in the view
of running a simulation pipeline, with realistic duration and number
of sources.

## Orbits

Orbits evolves along 2 dimensions: nSample x nLink with nLink << nSample.

Travel times are scalar, whereas positions are 3d vectors. 

The most efficient way of managing orbits is to make implicit loops
over time samples using ndarray. 

In the orbit object, time input/output argument is always considered
as a collection of time samples (not a single value).

Current software performances are: 

|                | Analytic    | 
| ---------------|-------------| 
| duration       | 1 yr        | 
| time step      | 5 s         | 
| position cpu   | 0.6 s       | 
| orbits memmax  | 800 Mb      |
| orbits memory  | 500 Mb      |
| nb spacecrafts | 3           |
| nb links       | 6           |

Alpha, cosine and sine of alpha are temporarily stored in memory in
order to save CPU time, as they are used 3 times when computing
positions. They are also reused for each link in projection.

## HpHc

HpHc evolves with source parameters, and time, that is nSource x nSample.

We use 2D ndarray to handle those 2 dimensions when computing
waveforms when possible (galatic binaries), to avoid explicit loop,
but this option is not used in the simulation pipeline.

The objects offers some method to merge or split HpHc object, such
that an HpHc object can contain one or several sources of the same
kind.

Below is a summary of assumptions and current software performance.


|               | MBHB        | EMRI  | GB   | SOBBH | 
| ------------- |-------------| ------|------|-------|
| duration      | 1 yr        | 1 yr  | 2 yr | 2 yr  |
| time step     | 5 s         | 15 s  | 15 s | 5 s   |
| hphc cpu      | 8 s         | 68 s  |0.3 s | 15 s  |
| interp        | 1.3s        | 0.5 s |1.9 s | 5.4 s |
| hphc memmax   | 1.4 Gb      |800 Mb |400 Mb| 1.2 Gb|
| hphc memory   | 400 Mb      |200 Mb |250 Mb| 600 Mb| 
| nb sources    | 100 ?       |10-100?|30e6  | 21e3  |


Numbers of given for a single source. 
MBHB memory consumption is high due to sampling in Fourier (5s).


## Projected strains

Projected strain dimensions of nSource x nLink x nSample. 

It uses quantities computed by Orbits which depends on time and link,
and HpHc wich depends on time and sources.

It is envisaged to dispatch the workload along number of sources, as
they can be treated independantly. The cost of doing this is a
duplication of orbits computation. A mitigated way would be to handle
batches of sources.

Therefore, the ProjectedStrain object uses implicit loop on time
samples, like orbits and hphc. 

Question left is: how to manage the loop over arms ? For now, this is
an explicit loop.

Below is a summary of assumptions and current software performance.

|                | MBHB        | EMRI  | GB   | SOBBH | 
| -------------- |-------------| ------|------|-------|
| duration       | 1 yr        | 1 yr  | 2 yr | 2 yr  |
| time step      | 5 s         | 15 s  | 15 s | 5 s   |
| proj cpu/comput| 100 s       | >5min |  8 s | 3 min |
| proj cpu/interp| 21 s        |  72 s | 11 s | 52 s  |
| proj memmax    | 2.25 Gb     |1.2 Gb |1.2 Gb| 4 Gb  | 
| nb sources     | 100 ?       |10-100?|30e6  | 21e3  |
| nb links       | 6           |  6    | 6    | 6     |


Numbers of given for a single source, projected on the 6 arms. 
Two options are studied:

- cpu/comput: hp,hc are computed twice for each link
- cpu/interp: hp,hc are computed once for all, on a slightly extended
  time range, and then interpolated twice for each link. The
  interpolator for hp and hc is computed once for all. 















