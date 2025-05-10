# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.2](https://github.com/Armavica/rebop/compare/v0.9.1...v0.9.2) - 2025-05-10

### Dependencies

- drop python 3.10 and upgrade dependencies
- upgrade rand to v0.9.1 and rand_distr to v0.5.1
- upgrade winnow to v0.7.10
- upgrade pyo3 to v0.24.2
- update rand to v0.9 and rand_distr to v0.5
- bump actions/create-github-app-token in the actions group

### Fixed

- allow inf or nan in parameter names (thank you @maurosilber!)
- add a few test cases for float parsing

### Miscellaneous

- update github CI

## [0.9.1](https://github.com/Armavica/rebop/compare/v0.9.0...v0.9.1) - 2025-01-23

### Added

- Allow parameters to be used in rate expressions.

### Changed

- change error message for species not involved in reactions
- allow species not in reactions to be init but with warning
- rename lib::Rate into lib::PRate
- rename PExpr::Concentration into PExpr::Variable

### Documentation

- fix parameters in docs index page
- add example with Michaelis-Menten rate
- improve docstrings to explain reaction rates
- fix some docstrings of the rust API

### Fixed

- use Mapping and Sequence instead of dict and list in type hints
- raise error on erroneous rate expressions
- make the rates more strict
- ensure that the dataset returned has a time coordinate

### Miscellaneous

- update pixi.lock with new version
- add test with crossed rates
- remove useless clone
- remove unused imports

## [0.9.0](https://github.com/Armavica/rebop/compare/v0.8.3...v0.9.0) - 2025-01-22

### Added

- Reaction rates can now be arbitrary formulas, instead of just LMA (Python
  API).
- Reaction rates can be internally sparse, leading to drastic speedups for
  models with many species (thank you @maurosilber!) (Python API).
- One can specify to not export all variables, with the `var_names` argument
  (thank you @maurosilber!) (Python API).

### Changed

- The `seed` argument was renamed to `rng` to conform to SPEC-7 (Python API).

### Dependencies

- upgrade pyo3 to 0.23

### Documentation

- organize badges on README
- add badges for the PyPI version and RTD docs
- fix name of optional dependencies in rtd config
- add rust tooling to the readthedocs environment
- fix readthedocs install
- add a readthedocs page
- document the Python interface

### Fixed

- fix arbitrary rates overwriting previous definitions
- fix benches with a sparsity argument

### Miscellaneous

- check all targets in pre-commit
- add tests for exprs
- robustify Gillespie tests
- Add tests for gillespie::Jump
- add tests for gillespie::Rate
- add tests for gillespie::Rate
- fix pixi build task
- update CI
- update pre-commit hooks
- upgrade dependencies
- cancel CI on repeated PRs
- enable dependabot update grouping
- use a markdown formatter
- add custom pre-commit CI message
- enable mypy
- update pre-commit hooks
- make ruff show fixes
- enable docstring format by ruff
- add pytest config options

## [0.8.3](https://github.com/Armavica/rebop/compare/v0.8.2...v0.8.3) - 2024-07-18

### Miscellaneous

- Use GitHub app to make release-plz trigger workflows

## \[0.8.2\] -- 2024-07-17

Test release with release-plz.

### Documentation

- Update image links in the README

## \[0.8.1\] -- 2024-07-16

Patch release to fix distribution to PyPI.

## \[0.8.0\] -- 2024-07-16

### Added

- Added a `advance_one_reaction` to the API, that allows calling the Python API
  `run` function with `nb_steps=0` in order to return all reactions.

### Changed

- Upgrade PyO3 to v0.22.1.

## \[0.7.0\] -- 2024-06-10

### Added

- It is now possible to specify a random seed from Python to run the simulator.

### Changed

- Upgraded PyO3 to v0.21.2.
- Dropped support for Python 3.9, according to the SPEC-0 schedule.

## \[0.6.1\] -- 2023-10-05

Patch release to fix distribution to PyPI. The release is now happening
automatically from GitHub CI.

## \[0.6.0\] -- 2023-10-05

### Added

- It is now possible to refer to species amounts in the macro rate expressions.
  The rate expressions are still multiplied by the law of mass action
  constructed from the reactions.
- Functions have been added to the API engine to set time and species amounts.
- An optional argument to the Python API has been added to specify reverse
  reactions more easily.

### Changed

- The Python results are now returned as an `xarray.Dataset` instead of a
  dictionary.
- The API engine has been rewritten to make it faster and more modular.
  Notably, the API to specify reaction rates in Rust changed from enums to
  vectors.

## \[0.5.1\] -- 2023-03-13

This release comes with the Python package `rebop` which exposes the Python
API.

### Changed

- The function `get_species` now takes a species instead of a reference to one.
- The libraries pyo3 and rand were upgraded to versions 0.15.1 and 0.8.5.

## \[0.5.0\] -- 2022-02-10

The major change of this release is a change in the macro DSL to write
reactions even more naturally, with `+` instead of `,`. Example:
`translation: mRNA => mRNA + GFP @ trans_rate`.

It comes also with some performance improvements, mostly on the API side (~5%).

### Added

- The Rngs can now be seeded.
- `SRate::LMA2` was reintroduced, along with `SRate::LMAn`, implementing
  2nd and nth-order LMA. They cannot be implemented in terms of
  `SRate::LMA`.
- Added nth-order reactants and products to the macro DSL.
- The libraries `rand` and `rand_distr` are now reexported to be more
  practical, notably this allows several system definitions in the
  same file.

### Changed

- Use `+` to separate species in macro DSL instead of `,`.

## \[0.4.0\] -- 2021-08-02

This release introduces a proof of concept for Python bindings, as well
as more performance improvements for the function-based API (20-35%).

### Added

- Python bindings.

### Changed

- `Gillespie` no longer has a const parameter. This breaking change
  was necessary for the definition of a system at run-time (for
  Python bindings).

## \[0.3.1\] -- 2021-08-01

Increased performance by 3-5% by optimizing exponential sampling.

## \[0.3.0\] -- 2021-08-01

Actual first release on [crates.io](https://crates.io).

### Changed

- Renamed the project from `bebop` to `rebop` because the name was taken.

## \[0.2.0\] -- 2021-08-01

First release on [crates.io](https://crates.io).

### Added

- Benchmark with other simulators.

### Removed

- `SRate::LMA2` which can already be constructed in terms of
  `SRate::LMA`.
- Removed the public visibility of `Rate::rate`, the function that
  evaluates the numerical value of the rate.

### Changed

- Renamed `choice` to `_choice` to make it clearer that it is an
  implementation macro.
- Relicensed `bebop` to MIT license.

## \[0.1.0\] -- 2021-07-31

First release.

### Added

- Two APIs for the simulation of chemical reaction networks:
  - `gillespie::Gillespie` which allows one to define a system at run
    time, and with several possible rate functions (Law of Mass Action
    but also Hill functions or Michaelis--Menten dynamics);
  - `define_system!` which requires that the system is defined at
    compilation. It only supports rate functions based on the Law
    of Mass Action. It runs typically three times as fast as the
    programmatic API.
