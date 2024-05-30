"""Class representing a standard LDC submission.

"""

import datetime
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any, Optional, List, Dict
import yaml

from ldc.eval.entry import LDCEntry, GBEntry


@dataclass
class LDCSubmission:
    """LDC submission data in full.

    Args:
        name (str): A name (ID) for the submission.
        entries (list[LDCEntry]): List of LDC entries contained in the submission.
        authors (list[str]): List of author names.
        emails (list[str]): List of author contact emails.
        challenge (str): LDC challenge that was attempted (e.g. "LDC1", "LDC2a").
        dataset (str): Specific LDC data filename that was analysed (e.g. "LDC1-3_VGB_v2.hdf5").
        date (datetime.date): Date when the submission was made.
        metadata (Optional[dict[str, Any]]): A dictionary with any additional metadata supplied.

    """

    name: str
    entries: List[LDCEntry] = field(repr=False)
    authors: List[str]
    emails: List[str]
    challenge: str
    dataset: str
    date: datetime.date
    metadata: Optional[Dict[str, Any]] = None

    @classmethod
    def from_yml(cls, path: Path):
        """Load an LDC submission that follows the standard format from a yaml file.

        For submissions that do not follow that standard, inherit from LDCSubmission
        and override this method

        Args:
            path (Path): Path to the .yml file that contains all submission information.

        Raises:
            NotImplementedError: This method cannot be implemented here until
        a standard LDC submission format is established.
        """
        raise NotImplementedError(
            "A standard LDC submission interface needs to be specified!"
        )

    def to_yml(self, path: Path) -> None:
        """Save an LDC submission to the standard yaml format, so it can be easily submitted.

        Args:
            path (Path): Path to the .yml file where the submission data will be saved.

        Raises:
            NotImplementedError: This method cannot be implemented until
        a standard LDC submission format is established.
        """

        raise NotImplementedError(
            "A standard LDC submission interface needs to be specified!"
        )

@dataclass
class SangriaGBSubmission(LDCSubmission):
    """ Sangria galactic binary submission.
    """
    @classmethod
    def from_yml(cls, path: Path):
        """Load an LDC submission that follows the standard format from a yaml file.

        Args:
            path (Path): Path to the .yml file that contains all submission information.
        """

        L = yaml.load(open(path, "r"), yaml.Loader)
        name = path
        authors = L["author"]
        df = pd.DataFrame(L["estimates"])

        # convert string into numbers, ignore symmetric error
        for i_n, n in enumerate(df.columns):
            v = [df[n][i].split(" +/- ") for i in range(len(df))]
            df[n] = [e[0] for e in v]
        entries = [GBEntry(i_n, df.loc[i_n]) for i_n in range(len(df))]
            
        emails = L["e-mail"]
        challenge = L["challenge"]
        dataset = L["dataset"]
        date = L["date"]
        return cls(name, entries, authors, emails, challenge, dataset, date)

