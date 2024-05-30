"""Class that contains the information in an LDC release.

This can be used to compare LDC submissions with "the truth".

"""
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

from ldc.common.series import TDI


@dataclass
class LDCDataset:
    """An LDC dataset.

    Should contain the catalog of "true" sources and TDI.

    Args:
        challenge (str): Name of the challenge (e.g. "LDC1", "LDC2a").
        version (int): Version of the dataset. Defaults to 1.
        tdi (Optional[TDI]): TDI data retrieved from the dataset.
        catalog (Optional[pd.DataFrame]): Dataframe with the true parameter values
        for each GW source in the dataset.
    """

    challenge: str
    version: int = 1
    tdi: Optional[TDI] = None
    catalog: Optional[pd.DataFrame] = None
    # TODO: how do we handle the Sangria (& future!) case where
    # there are multiple source types in the dataset?

    @property
    def name(self) -> str:
        """LISA data challenge name, version included.

        Returns:
            str: The name of this LDC.
        """
        return self.challenge + "_v" + str(self.version)

    @classmethod
    def from_hdf(cls, path: Path):
        """Load an LDC dataset from an hdf5 file.

        Args:
            path (Path): Path to the hdf5 file to the dataset.

        Raises:
            NotImplementedError: This method has not been implemented yet.
        """
        # TODO: Is the file structure of the dataset releases uniform between challenges?
        # Otherwise this class should be subclassed for each of Radler, Sangria, etc...
        raise NotImplementedError("Not implemented yet!")
