import xarray as xr
from .rebop import Gillespie, __version__

og_run = Gillespie.run

def run_xarray(self: Gillespie, init: dict[str, int], tmax: float, nb_steps: int) -> xr.Dataset:
    """
    Run the system until `tmax` with `nb_steps` steps.

    The initial configuration is specified in the dictionary `init`.
    Returns an xarray Dataset.
    """
    times, result = og_run(self, init, tmax, nb_steps)
    df = xr.Dataset(
        data_vars={
            name: xr.DataArray(values, dims="time", coords={"time": times})
            for name, values in result.items()
        }
    )
    return df

Gillespie.run = run_xarray
