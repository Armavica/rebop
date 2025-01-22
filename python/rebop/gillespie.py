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
        rate: float | str,
        reactants: list[str],
        products: list[str],
        reverse_rate: float | str | None = None,
    ) -> None:
        """Add a reaction to the system.

        Example
        -------
        ```python
        s = Gillespie()
        # Add the birth equation ø -> A at rate 14
        s.add_reaction(14, [], ["A"])
        # Add the reversible reaction A + B -> C
        # with forward rate 0.1 and reverse rate 0.01
        s.add_reaction(0.1, ["A", "B"], ["C"], 0.01)
        # Add the transformation B -> C with a non-LMA rate
        s.add_reaction("2.1 * B * C / (5 + C)", ["B"], ["C"])
        ```

        Parameters
        ----------
        rate : float | str
            If numeric, this is the rate constant of a Law of Mass Action.
            If string, this is the mathematical expression of the reaction rate.
        reactants : list[str]
            List of the species that are on the left-hand-side of the reaction.
        products : list[str]
            List of the species that are on the right-hand-side of the reaction.
        reverse_rate : float | str | None
            Rate of the reverse reaction, if this reaction is reversible.
            This is just a convenience to avoid using `add_reaction` a second
            time for the reverse reaction.
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
        sparse: bool | None = None,
        var_names: Sequence[str] | None = None,
    ) -> xr.Dataset:
        """Run the system until `tmax` with `nb_steps` steps.

        Parameters
        ----------
        init : dict[str, int]
            Dictionary that indicates the initial condition. It can omit
            species, their initial number will be 0.
        tmax : float
            Simulation end time.
        nb_steps : int
            Number of steps to return, equally distributed between 0 and `tmax`.
            If 0, then all reactions are returned.
        rng : RNGLike | SeedLike | None
            Numpy `Generator`, `BitGenerator` or seed.
        sparse : bool | None
            Whether to internally represent reactions as dense or sparse.
            Sparse reactions will tend to be faster for systems with many
            species. If `None`, rebop will choose depending on heuristics.
        var_names : Sequence[str] | None
            List of variable names to return, if one does not want to save
            all variables. If `None`, all variables will be saved.

        Returns
        -------
        xr.Dataset
            `xarray.Dataset` object with a `time` dimension and variables
            being chemical species (all of them, or just the subset defined
            by `var_names`).
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
