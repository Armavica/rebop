# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
    * `gillespie::Gillespie` which allows a system to be defined at run
        time, and with several possible rate functions (Law of Mass Action
        but also Hill functions or Michaelis--Menten dynamics);
    * `define_system!` which requires that the system is defined at
        compilation.  It only supports rate functions based on the Law
        of Mass Action.  It runs typically three times as fast as the
        programmatic API.