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
def test_sir(seed: int | None) -> None:
    sir = sir_model()
    ds = sir.run({"S": 999, "I": 1}, tmax=250, nb_steps=250, seed=seed)
    assert isinstance(ds, xr.Dataset)
    npt.assert_array_equal(ds.time, np.arange(251))
    assert all(ds.S >= 0)
    assert all(ds.I >= 0)
    assert all(ds.R >= 0)
    assert all(ds.S <= 999)
    assert all(ds.I <= 1000)
    assert all(ds.R <= 1000)
    npt.assert_array_equal(ds.S + ds.I + ds.R, [1000] * 251)


def test_fixed_seed() -> None:
    sir = sir_model()
    ds = sir.run({"S": 999, "I": 1}, tmax=250, nb_steps=250, seed=42)

    assert ds.S[-1] == 0
    assert ds.I[-1] == 166
    assert ds.R[-1] == 834


@pytest.mark.parametrize("seed", range(10))
def test_all_reactions(seed: int) -> None:
    tmax = 250
    sir = sir_model()
    ds = sir.run({"S": 999, "I": 1}, tmax=tmax, nb_steps=0, seed=seed)
    assert ds.time[0] == 0
    assert ds.time[-1] > tmax
    assert all(ds.time.diff(dim="time") > 0)
    if np.isinf(ds.time[-1].to_numpy()):
        ds = ds.isel(time=range(ds.time.size - 1))
    dds = ds.diff(dim="time")
    assert set(dds.S.to_numpy()) <= {-1, 0}
    assert set(dds.I.to_numpy()) <= {-1, 1}
    assert set(dds.R.to_numpy()) <= {0, 1}


def test_dense_vs_sparse() -> None:
    sir = sir_model()
    init = {"S": 999, "I": 1}
    tmax = 250
    nb_steps = 250
    seed = 42
    ds_dense = sir.run(init, tmax=tmax, nb_steps=nb_steps, seed=seed, sparse=False)
    ds_sparse = sir.run(init, tmax=tmax, nb_steps=nb_steps, seed=seed, sparse=True)
    assert (ds_dense == ds_sparse).all()


@pytest.mark.parametrize("nb_steps", [0, 250])
def test_var_names(nb_steps: int) -> None:
    all_variables = {"S", "I", "R"}
    subset_to_save = ["S", "I"]
    remaining = all_variables.difference(subset_to_save)

    sir = sir_model()
    init = {"S": 999, "I": 1}
    tmax = 250
    seed = 0

    ds_all = sir.run(init, tmax=tmax, nb_steps=nb_steps, seed=seed, var_names=None)
    ds_subset = sir.run(
        init, tmax=tmax, nb_steps=nb_steps, seed=seed, var_names=subset_to_save
    )

    for s in subset_to_save:
        assert s in ds_subset
    for s in remaining:
        assert s not in ds_subset

    assert ds_all[subset_to_save] == ds_subset
