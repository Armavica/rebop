from __future__ import annotations

import xarray as xr

from .rebop import Gillespie, __version__  # type: ignore[attr-defined]

__all__ = ("Gillespie", "__version__")

og_run = Gillespie.run


def run_xarray(  # noqa: PLR0913 too many parameters in function definition
    self: Gillespie,
    init: dict[str, int],
    tmax: float,
    nb_steps: int,
    seed: int | None = None,
    *,
    sparse: bool = False,
) -> xr.Dataset:
    """Run the system until `tmax` with `nb_steps` steps.

    The initial configuration is specified in the dictionary `init`.
    Returns an xarray Dataset.
    """
    times, result = og_run(self, init, tmax, nb_steps, seed, sparse=sparse)
    ds = xr.Dataset(
        data_vars={
            name: xr.DataArray(values, dims="time", coords={"time": times})
            for name, values in result.items()
        },
    )
    return ds


Gillespie.run = run_xarray  # type: ignore[method-assign]
