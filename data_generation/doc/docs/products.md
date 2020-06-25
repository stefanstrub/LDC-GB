# The pipeline data products

## GW source data products

For each kind of source `{source}`, the pipeline produces:

- a catalog of selected source: `{source}.npy`
- a projected strain file, which contains the sum of all sources of this kind: `{source}-y.h5`

The merged projected strain file is named `sum-y.h5`. 

## L1 data products

L0 data are not saved on disk for now, due to issues with the HDF5
writer node in LISANode (to do so, we would need to store all the data
in memory first, which is not possible for long duration run). 

L1 data are produced by LISANode, for each type of source (without
instrumental noise), and for the sum of them (with instrumental noise).
They are saved in the `{source}-lisanode/{source-tdi.h5}` files. 
