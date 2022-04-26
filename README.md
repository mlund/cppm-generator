[![Rust](https://github.com/mlund/cppm-generator/actions/workflows/rust.yml/badge.svg)](https://github.com/mlund/cppm-generator/actions/workflows/rust.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6487965.svg)](https://doi.org/10.5281/zenodo.6487965)

# CPPM generator

## Overview

This is a single command-line tool that generates spherical
Charged Patchy Particle Models (CPPM) which mimic bio-colloidal
particles, _e.g._ globular proteins.
The figure below shows examples of interacting CPPMs taken from https://doi.org/10.1063/1.4928077:

<img src="https://aip.scitation.org/action/showOpenGraphArticleImage?doi=10.1063/1.4928077&id=images/medium/1.4928077.figures.f14.gif" height="100" />

`cppm-generator` generates CPPMs by placing neutral, positive, and negative particles
on the surface of a sphere, and minimise the (free) energy using
Metropolis-Hastings Monte Carlo sampling.

## Command line usage:

The default parameters produce an isotropic, charged particle similar to `P00` from Table 1 in
[this](https://doi.org/10.48550/arXiv.1701.02457) publication.
It is also possible to impose a target molecular dipole moment using the `--dipole` option.

~~~
$ cppm-generator --help

USAGE:
    cppm-generator [OPTIONS] --file <FILE>

OPTIONS:
    -b, --bjerrum-length <BJERRUM_LENGTH>    Bjerrum length (Å) [default: 7.0]
    -h, --help                               Print help information
    -m, --minus <NUM_MINUS>                  Number of negative (-1e) particles [default: 37]
    -N <NUM_TOTAL>                           Total number of particles [default: 643]
    -o, --file <FILE>                        Output structure (.xyz or .pqr)
    -p, --plus <NUM_PLUS>                    Number of positive (+1e) particles [default: 29]
    -r, --radius <RADIUS>                    Sphere radius (Å) [default: 20.0]
    -s, --steps <STEPS>                      Number of Monte Carlo iterations [default: 10000]
    -u, --dipole <TARGET_DIPOLE_MOMENT>      Target dipole moment (Debye)
    -V, --version                            Print version information
~~~

## Details and status:

- [x] Random walk on a sphere using spherical coordinates
- [x] Particle-particle interactions using a Coulomb/softcore potential
- [x] Arbitrary mixing of neutral and charged particles
- [x] Output to `.xyz` and `.pqr` files
- [x] Command line interface
- [x] Dipole moment analysis
- [ ] External electric field to induce arbitrary patches
- [x] Constrain to target dipole moment w. harmonic potential
- [x] Written in Rust
- [ ] Use [uon](https://crates.io/crates/uom) for dimensional analysis
- [x] IO error handling
- [ ] Unittests
- [ ] Logging support

## Building and running

Building from source requires a [Rust installation](https://www.rust-lang.org/tools/install).

~~~ bash
cd cppm-generator/
cargo build --release
./target/release/cppm-generator --help

# Optionally install in $HOME/.cargo/bin:
cargo install --path .
~~~

## Motivation

Besides some scientific use, this project is mainly a first dive
into the Rust programming language.
