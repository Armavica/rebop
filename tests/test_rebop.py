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
    ds = sir.run({"S": 999, "I": 1}, tmax=250, nb_steps=250, rng=seed)
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
    ds = sir.run({"S": 999, "I": 1}, tmax=250, nb_steps=250, rng=42)

    assert ds.S[-1] == 0
    assert ds.I[-1] == 227
    assert ds.R[-1] == 773


@pytest.mark.parametrize("seed", range(10))
def test_all_reactions(seed: int) -> None:
    tmax = 250
    sir = sir_model()
    ds = sir.run({"S": 999, "I": 1}, tmax=tmax, nb_steps=0, rng=seed)
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
    ds_auto = sir.run(init, tmax=tmax, nb_steps=nb_steps, rng=seed, sparse=None)
    ds_dense = sir.run(init, tmax=tmax, nb_steps=nb_steps, rng=seed, sparse=False)
    ds_sparse = sir.run(init, tmax=tmax, nb_steps=nb_steps, rng=seed, sparse=True)
    assert (ds_dense == ds_auto).all()
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

    ds_all = sir.run(init, tmax=tmax, nb_steps=nb_steps, rng=seed, var_names=None)
    ds_subset = sir.run(
        init, tmax=tmax, nb_steps=nb_steps, rng=seed, var_names=subset_to_save
    )

    for s in subset_to_save:
        assert s in ds_subset
    for s in remaining:
        assert s not in ds_subset

    assert ds_all[subset_to_save] == ds_subset


def test_arbitrary_rates() -> None:
    s = rebop.Gillespie()
    with pytest.raises(ValueError, match="Rate expression not understood"):
        s.add_reaction("+", [], ["A"])
    with pytest.raises(ValueError, match="Rate expression not understood"):
        s.add_reaction("1+", [], ["A"])

    s = rebop.Gillespie()
    s.add_reaction("B", [], ["A"])
    with pytest.warns(
        UserWarning,
        match="species are not involved in any reaction",
    ):
        ds = s.run({"B": 1}, tmax=10, nb_steps=100)
    assert ds.A[-1] >= 1
    np.testing.assert_equal(ds.B.values, [1] * 101)

    s = rebop.Gillespie()
    s.add_reaction("B", [], ["A"])
    s.add_reaction(1, [], ["B"])
    ds_with_b = s.run({}, tmax=10, nb_steps=100)

    assert ds_with_b.A[-1] >= 1


def test_arbitrary_rates_crossed() -> None:
    s = rebop.Gillespie()
    s.add_reaction("B", [], ["A"])
    s.add_reaction("A", [], ["B"])
    ds = s.run({}, tmax=10, nb_steps=10)
    expected = xr.Dataset(
        {"A": ("time", [0] * 11), "B": ("time", [0] * 11)},
        coords={"time": np.linspace(0, 10, 11)},
    )
    xr.testing.assert_identical(ds, expected)

    ds = s.run({"A": 1}, tmax=10, nb_steps=10)
    assert ds.A[-1] > 1
    assert ds.B[-1] > 0

    ds = s.run({"B": 1}, tmax=10, nb_steps=10)
    assert ds.A[-1] > 0
    assert ds.B[-1] > 1


def test_arbitrary_rates_2() -> None:
    s = rebop.Gillespie()
    s.add_reaction(14, [], ["A"])
    s.add_reaction(0.1, ["A", "B"], ["C"], 0.01)
    s.add_reaction("0.2 * B * C / (5 + C)", ["B"], ["D"])
    ds = s.run({"B": 1000}, tmax=100, nb_steps=100)

    assert (ds.B + ds.C + ds.D == 1000).all()
    assert ds.D[-1] >= 1


def test_run_empty() -> None:
    s = rebop.Gillespie()
    ds = s.run({}, tmax=10, nb_steps=10)
    expected = xr.Dataset(coords={"time": np.linspace(0, 10, 11)})
    xr.testing.assert_identical(ds, expected)


def test_parameters() -> None:
    s = rebop.Gillespie()
    s.add_reaction(4, ["A"], ["B"])
    with pytest.raises(ValueError, match="Species B cannot also be a parameter"):
        s.run({}, 10, 10, params={"B": 4.2})

    s = rebop.Gillespie()
    s.add_reaction("k", [], ["A"])
    with pytest.raises(ValueError, match="Parameter k should have a value"):
        s.run({}, 10, 10)
    ds = s.run({}, 10, 10, params={"k": 0.4})
    assert ds.A[-1] > 0
