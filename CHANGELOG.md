# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.6.1] -- 2023-10-05

Patch release to fix distribution to PyPI.  The release is now happening
automatically from GitHub CI.

## [0.6.0] -- 2023-10-05

### Added

* It is now possible to refer to species amounts in the macro rate expressions.
  The rate expressions are still multiplied by the law of mass action
  constructed from the reactions.
* Functions have been added to the API engine to set time and species amounts.
* An optional argument to the Python API has been added to specify reverse
  reactions more easily.

### Changed

* The Python results are now returned as an `xarray.Dataset` instead of a
  dictionary.
* The API engine has been rewritten to make it faster and more modular.
  Notably, the API to specify reaction rates in Rust changed from enums to
  vectors.

## [0.5.1] -- 2023-03-13

This release comes with the Python package `rebop` which exposes the Python
API.

### Changed

* The function `get_species` now takes a species instead of a reference to one.
* The libraries pyo3 and rand were upgraded to versions 0.15.1 and 0.8.5.

## [0.5.0] -- 2022-02-10

The major change of this release is a change in the macro DSL to write
reactions even more naturally, with `+` instead of `,`.  Example:
`translation: mRNA => mRNA + GFP @ trans_rate`.

It comes also with some performance improvements, mostly on the API side (~5%).

### Added

* The Rngs can now be seeded.
* `SRate::LMA2` was reintroduced, along with `SRate::LMAn`, implementing
    2nd and nth-order LMA.  They cannot be implemented in terms of
    `SRate::LMA`.
* Added nth-order reactants and products to the macro DSL.
* The libraries `rand` and `rand_distr` are now reexported to be more
    practical, notably this allows several system definitions in the
    same file.

### Changed

* Use `+` to separate species in macro DSL instead of `,`.

## [0.4.0] -- 2021-08-02

This release introduces a proof of concept for Python bindings, as well
as more performance improvements for the function-based API (20-35%).

### Added

* Python bindings.

### Changed

* `Gillespie` no longer has a const parameter.  This breaking change
    was necessary for the definition of a system at run-time (for
    Python bindings).

## [0.3.1] -- 2021-08-01

Increased performance by 3-5% by optimizing exponential sampling.

## [0.3.0] -- 2021-08-01

Actual first release on [crates.io](https://crates.io).

### Changed

* Renamed the project from `bebop` to `rebop` because the name was taken.

## [0.2.0] -- 2021-08-01

First release on [crates.io](https://crates.io).

### Added

* Benchmark with other simulators.

### Removed

* `SRate::LMA2` which can already be constructed in terms of
    `SRate::LMA`.
* Removed the public visibility of `Rate::rate`, the function that
    evaluates the numerical value of the rate.

### Changed

* Renamed `choice` to `_choice` to make it clearer that it is an
    implementation macro.
* Relicensed `bebop` to MIT license.

## [0.1.0] -- 2021-07-31

First release.

### Added

* Two APIs for the simulation of chemical reaction networks:
    * `gillespie::Gillespie` which allows one to define a system at run
        time, and with several possible rate functions (Law of Mass Action
        but also Hill functions or Michaelis--Menten dynamics);
    * `define_system!` which requires that the system is defined at
        compilation.  It only supports rate functions based on the Law
        of Mass Action.  It runs typically three times as fast as the
        programmatic API.
