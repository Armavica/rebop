"""rebop is a fast stochastic simulator for well-mixed chemical reaction networks."""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
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
        reactants: Sequence[str],
        products: Sequence[str],
        reverse_rate: float | str | None = None,
    ) -> None:
        """Add a reaction to the system.

        Reaction rates can be specified with a number (for a reaction obeying
        the Law of Mass Action) or a string (for an arbitrary reaction rate).
        The string can involve any species involved in reactions, and also
        parameters defined with the `params` argument of the `run` method.

        If you can, use the law of mass action, which is probably going to be
        more efficient and can be more correct. For example, for a dimerisation
        equation:

        ```python
        s = Gillespie()
        # Correct, and recommended
        s.add_reaction(4.2, ["A", "A"], ["AA"])
        # Correct, but not recommended
        s.add_reaction("4.2 * A * (A-1)", ["A", "A"], ["AA"])
        # Incorrect (this would be a constant propensity)
        s.add_reaction("4.2", ["A", "A"], ["AA"])
        # Incorrect (incorrect expression)
        s.add_reaction("4.2 * A^2", ["A", "A"], ["AA"])
        ```

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
        s.add_reaction("2.1 * B * A / (5 + A)", ["B"], ["C"])
        # Add the transformation B -> C with a non-LMA rate and parameters
        # that need to be defined by passing `params={"k": 2.1, "K_A": 5}` to
        # the method `run`
        s.add_reaction("k * B * A / (K_A + A)", ["B"], ["C"])
        ```

        Parameters
        ----------
        rate : float | str
            If numeric, this is the rate constant of a Law of Mass Action.
            If string, this is the mathematical expression of the reaction rate.
        reactants : Sequence[str]
            List of the species that are on the left-hand-side of the reaction.
        products : Sequence[str]
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
        init: Mapping[str, int],
        tmax: float,
        nb_steps: int,
        *,
        params: Mapping[str, float] | None = None,
        rng: RNGLike | SeedLike | None = None,
        sparse: bool | None = None,
        var_names: Sequence[str] | None = None,
    ) -> xr.Dataset:
        """Run the system until `tmax` with `nb_steps` steps.

        Parameters
        ----------
        init : Mapping[str, int]
            Dictionary that indicates the initial condition. It can omit
            species, their initial number will be 0.
        tmax : float
            Simulation end time.
        nb_steps : int
            Number of steps to return, equally distributed between 0 and `tmax`.
            If 0, then all reactions are returned.
        params : Mapping[str, float] | None
            Dictionary of values for the parameters that appear in the rates.
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
        params = {} if params is None else params
        rng_ = np.random.default_rng(rng)
        seed = rng_.integers(np.iinfo(np.uint64).max, dtype=np.uint64)
        try:
            self.gillespie.set_init(init)
        except UserWarning as e:
            warnings.warn(e, stacklevel=2)
        times, result = self.gillespie.run(
            tmax,
            nb_steps,
            params=params,
            seed=seed,
            sparse=sparse,
            var_names=var_names,
        )
        ds = xr.Dataset(
            data_vars={
                name: xr.DataArray(values, dims="time", coords={"time": times})
                for name, values in result.items()
            },
            coords={"time": times},
        )
        return ds
