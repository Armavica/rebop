from __future__ import annotations

from collections.abc import Sequence
from typing import TypeAlias

import numpy as np
import xarray as xr

from .rebop import Gillespie, __version__  # type: ignore[attr-defined]

SeedLike: TypeAlias = int | np.integer | Sequence[int] | np.random.SeedSequence
RNGLike: TypeAlias = np.random.Generator | np.random.BitGenerator


__all__ = ("Gillespie", "__version__")


def run_xarray(  # noqa: PLR0913 too many parameters in function definition
    self: Gillespie,
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
    times, result = self._run(
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


Gillespie.run = run_xarray  # type: ignore[method-assign]
