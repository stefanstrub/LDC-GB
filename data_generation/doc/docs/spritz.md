# Spritz production

The spritz data set consists in 3 data files:

- mbhb1: a loud event with 3 short glitches
- mbhb2: a long glitch on top of merger
- vgb: verification binairies with LPF-like glitch distribution

The data contains:

- glitches
- gaps
- non stationary noise which is the background galactic binaries signal (ie unresolved GB)

<img src="spritz/dag.png"  alt="Pipeline DAG" width="750"/> 

## Ancillary data

To run the spritz pipeline, one needs the following data files:

- `VGB.h5` for the 17 verification binaries;
- `Q3d_complete` for the 692 MBHB;
- `2021-04-15-intervals_ordinary.txt` which set the time interval between glitches for VGB data set
- `checkpoint_lr_90001_512.pt` which set the parameter distribution of glitches for VGB data set
- `lpf-glitch-library.h5` the LPF glitch template library
- `reduced_1_yr.h5` the reduced galaxy computed by FoM pipeline over 1 year of observation

