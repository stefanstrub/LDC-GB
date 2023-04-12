LDC.2a Sangria

Root directory is /lisa/sangria

- dgb-pipeline/run1, igb-pipeline/run1, run2: data released on the 2020/10/08 see https://lisa-ldc.lal.in2p3.fr/file
- dgb-pipeline/run1, igb-pipeline/run1, run3: change interpolation of dgb and igb with order=3->9, fix fdot in 3 vgb and ampl in 1. 
- dgb-pipeline/run1, igb-pipeline/run1, run4: as run3 but change interpolation of dgb and igb with 0-padding in Fourier domain
- dgb-pipeline/run2, igb-pipeline/run1, run5: as run3 + new dgb catalog with SNR cut and uniform inclination distribution
- dgb-pipeline/run2, igb-pipeline/run1, run6: as run4 + new dgb catalog with SNR cut and uniform inclination distribution
- dgb-pipeline/run2, igb-pipeline/run1, v2: as run6 + new vgb catalog and mbhb sorted by tc

Note:

Since oct-20 release, we try to address several issues reported by US team:
- distribution of inclination in dgb catalog
- artifact at high freq in y_t and tdi (above max gb freq of 0.04 Hz) due to dgb upsampling
- very high amplitude dgb source (snr > 5000), which is making the above artifact worst

We also change the pipeline to follow latest lisanode developments and have coherent output between sangria and spritz:
- orbits from file + upsampling
- arm projection use orbits from file
- input gw strain at dt=2.5s (was dt=3s before)
- lisanode at 4Hz (was 3Hz before)
- tdi downsampling out of lisanode, to reduce compilation time
- tdi downsampling use Kaiser filter (was Elliptic before)
- output file contains noise free tdi without negative time samples

