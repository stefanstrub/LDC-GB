"""Different LDC entry types.

"""
from abc import ABC
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, List

import corner
import matplotlib.pyplot as plt
import pandas as pd

from ldc.common.series import TDI

# Type definition: A function that takes in a pandas Series (parameters)
# and returns a TDI object, is called a Simulator.
Simulator = Callable[[pd.Series], TDI]


@dataclass
class LDCEntry(ABC):
    """Abstract type representing an LDC entry of any kind.

    Any methods defined in this class will work for all LDC entry types, unless overloaded.
    Specific LDC entry types (e.g. MBHB, GB) should inherit from this class
    and implement/replace what needs to be specific for them.

    Args:
        name (str): A name (ID) for this specific entry.
        point_estimate (Optional[pd.Series]): Parameter values of the obtained point estimate.
        samples_path (Optional[Path]): Path to the file containing posterior samples for this entry.
        load_data (Callable[[Path], pd.DataFrame]): Function that reads the posterior file
        and returns the samples. Defaults to pd.read_csv without arguments.
        simulate (Optional[Simulator]): function that takes in the parameters of the entry
        and returns a TDI object with simulator data.
        cache_samples (bool): If True, load samples and keep them in memory. Defaults to False.
        cache_tdi (bool): If True and both point_estimate and simulate are provided,
        precompute TDI of the point estimate and keep it in memory. Defaults to False.

    """

    name: str
    point_estimate: Optional[pd.Series] = None
    samples_path: Optional[Path] = None
    load_data: Callable[[Path], pd.DataFrame] = pd.read_csv
    simulate: Optional[Simulator] = None
    cache_samples: bool = field(default=False, repr=False)
    cache_tdi: bool = field(default=False, repr=False)
    _cached_samples: Optional[pd.DataFrame] = field(
        default=None, init=False, repr=False
    )
    _cached_tdi: Optional[TDI] = field(default=None, init=False, repr=False)
    # TODO: what to do when submissions provide one- or two-sided uncertainties
    # with their point estimates?

    def __post_init__(self):
        if self.cache_samples:
            self._cached_samples = self.load_data(self.samples_path)
            if self.cache_tdi and self.point_estimate and self.simulate:
                self._cached_tdi = self.simulate(self.point_estimate)

    @property
    def samples(self) -> pd.DataFrame:
        """Retrieve posterior samples.

        Returns:
            pd.DataFrame: Posterior samples corresponding to this entry.
        """
        if self._cached_samples:
            return self._cached_samples
        return self.load_data(self.samples_path)

    @property
    def tdi(self) -> TDI:
        """Get TDI of the provided point estimate for the entry.

        Returns:
            TDI: TDI corresponding to this entry.
        """
        if self._cached_tdi:
            return self._cached_tdi
        return self.simulate(self.point_estimate)

    def plot_corner(
        self,
        columns: Optional[List[str]] = None,
        fig: Optional[plt.Figure] = None,
        **corner_kwargs,
    ) -> plt.Figure:
        """Perform a basic corner plot.

        Args:
            columns (Optional[list[str]], optional): List of parameters to include in the plot.
            If None, plot all. Defaults to None.
            fig (Optional[plt.Figure], optional): Matplotlib Figure in which to perform
            the corner plot. Defaults to None.
            **corner_kwargs: Additional parameters that will be passed to corner.corner()
        Returns:
            plt.Figure: _description_
        """
        if fig is None:
            fig = plt.gcf()

        samples = self.samples
        if columns:
            samples = samples[columns]

        corner_kwargs = {
            "bins": 40,
            "quantiles": [0.05, 0.5, 0.95],
            "labels": samples.columns,
            "hist_kwargs": {"density": True, "lw": 2},
            "plot_datapoints": False,
            "fill_contours": False,
            "show_titles": False,
            "title_fmt": ".3g",
            "use_math_test": True,
            **corner_kwargs,
        }

        fig = corner.corner(samples, fig=fig, **corner_kwargs)
        return fig


@dataclass
class MBHBEntry(LDCEntry):
    """LDC entry representing a found MBHB source.

    All tools/visualizations specific to MBHB sources should be defined as methods in this class.

    """


@dataclass
class GBEntry(LDCEntry):
    """LDC entry representing a found GB source.

    All tools/visualizations specific to GB sources should be defined as methods in this class.

    """

    def __post_init__(self):
        pass


    
