# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
