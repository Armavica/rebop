import numpy as np
import numpy.testing as npt
import pytest
import rebop
import xarray as xr


def sir_model(transmission: float = 1e-4, recovery: float = 0.01) -> rebop.Gillespie:
    sir = rebop.Gillespie()
    sir.add_reaction(transmission, ["S", "I"], ["I", "I"])
    sir.add_reaction(recovery, ["I"], ["R"])
    return sir


@pytest.mark.parametrize("seed", [None, *range(10)])
def test_sir(seed: int):
    sir = sir_model()
    ds = sir.run({"S": 999, "I": 1}, tmax=250, nb_steps=250, seed=seed)
    assert isinstance(ds, xr.Dataset)
    npt.assert_array_equal(ds.time, np.arange(251))
    assert all(ds.S >= 0)
    assert all(ds.I >= 0)
    assert all(ds.R >= 0)
    assert all(ds.S <= 1000)
    assert all(ds.I <= 1000)
    assert all(ds.R <= 1000)
    npt.assert_array_equal(ds.S + ds.I + ds.R, [1000] * 251)


def test_fixed_seed():
    sir = sir_model()
    ds = sir.run({"S": 999, "I": 1}, tmax=250, nb_steps=250, seed=42)

    assert ds.S[-1] == 0
    assert ds.I[-1] == 166
    assert ds.R[-1] == 834
