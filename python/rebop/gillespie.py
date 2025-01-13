"""rebop is a fast stochastic simulator for well-mixed chemical reaction networks."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TypeAlias

import numpy as np
import xarray as xr

from rebop import _lib

__all__ = ("Gillespie", "__version__")


__version__: str = _lib.__version__

SeedLike: TypeAlias = int | np.integer | Sequence[int] | np.random.SeedSequence
RNGLike: TypeAlias = np.random.Generator | np.random.BitGenerator


class Gillespie:
    """Reaction system composed of species and reactions."""

    def __init__(self) -> None:
        """Initialize a solver."""
        self.gillespie = _lib.Gillespie()

    def add_reaction(
        self,
        rate: float,
        reactants: list[str],
        products: list[str],
        reverse_rate: float | None = None,
    ) -> None:
        """Add a Law of Mass Action reaction to the system.

        The forward reaction rate is `rate`, while `reactants` and `products`
        are lists of respectively reactant names and product names.
        Add the reverse reaction with the rate `reverse_rate` if it is not
        `None`.
        """
        self.gillespie.add_reaction(rate, reactants, products, reverse_rate)

    def __str__(self) -> str:
        """Return a textual representation of the reaction system."""
        return str(self.gillespie)

    def run(  # noqa: PLR0913 too many parameters in function definition
        self,
        init: dict[str, int],
        tmax: float,
        nb_steps: int,
        *,
        rng: RNGLike | SeedLike | None = None,
        sparse: bool = False,
        var_names: Sequence[str] | None = None,
    ) -> xr.Dataset:
        """Run the system until `tmax` with `nb_steps` steps.

        The initial configuration is specified in the dictionary `init`.
        Returns an xarray Dataset.
        """
        rng_ = np.random.default_rng(rng)
        seed = rng_.integers(np.iinfo(np.uint64).max, dtype=np.uint64)
        times, result = self.gillespie.run(
            init,
            tmax,
            nb_steps,
            seed=seed,
            sparse=sparse,
            var_names=var_names,
        )
        ds = xr.Dataset(
            data_vars={
                name: xr.DataArray(values, dims="time", coords={"time": times})
                for name, values in result.items()
            },
        )
        return ds
