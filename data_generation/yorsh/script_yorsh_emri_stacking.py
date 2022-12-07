""" This script stacks several user-specified fields from a single .hdf5 file
and store the resulting stacked TDIs XYZ in a new .hdf5 file. It is first
dedicated to the Yorsh dataset and more specically the FastSchwarschildFluxEccentric
model for EMRIs.
"""

import argparse
import h5py as h5
import numpy as np

# DICT_FIELDS = [
#     'emri1': None,
#     'emri2': None,
#     'emri3': None,
#     'inst_noise': None,
#     'gal_noise': None,
# ]

# # NOTE: remains to be completed then tested !
# def stack(dir, fields):
#     N = 1
#     Ntemp = 1
#     # read data
#     with h5.File(XXX) as f_:
#         for field in DICT_FIELDS:
#             DICT_FIELDS[field] = np.vstack(
#                 (
#                     f_.get(f'{field}/t'), 
#                     f_.get(f'{field}/X'), 
#                     f_.get(f'{field}/Y'), 
#                     f_.get(f'{field}/Z')
#                 )
#             ) # array of size (N, 4)
#             if (N==1) or (Ntemp < N):
#                 N = DICT_FIELDS[field].shape[0]
#                 dtype = DICT_FIELDS[field].dtype
#     # stack non None keys in DICT_FIELDS
#     stacked_data = np.zeros((N, 4), dtype=dtype)
#     for field in DICT_FIELDS:
#         if DICT_FIELDS[field] is not None:
#             stacked_data += DICT_FIELDS[field]
#     return stacked_data


def stack(fn: str, fields: list) -> dict:
    """Stack TDIs XYZ from various specified fields.

    Args:
        fn (string): path to input .hdf5 file containing TDIs XYZ to be merged.
        Expected structure of the .hdf5 file is described blow. Note that a 't'
        field is required.

            input.hdf5
               |- emri1
                  |- X
                  |- Y
                  |- Z
               |- emri2
                  |- X
                  |- Y
                  |- Z
               |- ...
               |- t

        fields (list): list of fields to be merged. Expected structure 
        (for example) is:
        
            ['emri1', 'emri2', 't'].

    Raises:
        ValueError: not enough fields have been provided.

    Returns:
        dict: dictionnary containing stacked TDIs. Has the following structure:
        {'t': ..., 'X': ..., 'Y': ..., 'Z': ...} where X, Y and Z are the 
        stacked TDIs. 
    """
    # sanity checks on inputs
    # At least one field is required in plus of the 
    # already expected time field...
    if not len(fields) > 1:
        raise ValueError(f"At least one field (in plusof time field) is required (found {len(fields)})")
    # open .hdf5 file
    with h5.File(fn) as f_:
        # read extensions names
        # that are not time
        extnames = [ e for e in list(f_.keys()) if e != 't' ]
        # read and save time vector
        t = f_['t']
        # declare final dictionnary
        stacked = {
            't': t,
            'X': np.zeros(t.size, dtype=np.float64),
            'Y': np.zeros(t.size, dtype=np.float64),
            'Z': np.zeros(t.size, dtype=np.float64),
        }
        # for all fields in input arg
        # stack X, Y, Z subfields
        for ext in extnames:
            if ext in fields:
                stacked['X'] += f_[ext]['X'][:]
                stacked['Y'] += f_[ext]['Y'][:]
                stacked['Z'] += f_[ext]['Z'][:]
    return stacked


if __name__ == "__main__":
    # user interface
    args = argparse.ArgumentParser()
    args.add_argument('-i', '--input', type=str, help="Path to where input .hdf5 file is stored")
    args.add_argument('-f', '--fields', type=list, help="List of desired fields to be stacked")

    # stack data extracted from required fields
    stacked_data = stack(args.i, args.f)

    # write
    name = 'stacked.h5'
    with h5.File(name, 'w') as f_:
        g = f_.create_group(name.lower())
        g.create_dataset('t', data=stacked_data["t"])
        g.create_dataset('Xstacked', data=stacked_data["X"])
        g.create_dataset('Ystacked', data=stacked_data["Y"])
        g.create_dataset('Zstacked', data=stacked_data["Z"])


