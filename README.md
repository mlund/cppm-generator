[![Rust](https://github.com/mlund/cppm-generator/actions/workflows/rust.yml/badge.svg)](https://github.com/mlund/cppm-generator/actions/workflows/rust.yml)

# `cppm-generator`

## Overview

This is a small command-line tool to generate spherical
_Charged Patchy Particle Models_ (CPPM) which mimics bio-colloidal
particles, _e.g._ globular proteins.
The figure below shows examples of interacting CPPMs taken from https://doi.org/10.1063/1.4928077:

<img src="https://aip.scitation.org/action/showOpenGraphArticleImage?doi=10.1063/1.4928077&id=images/medium/1.4928077.figures.f14.gif" height="100" />

This tool generates CPPMs by placing neutral, positive, and negative particles
on the surface of a sphere, and minimise the (free) energy using
Metropolis-Hastings Monte Carlo sampling.

## Command line usage:

The default parameters produces an isotropic, charged particle simular to `P00` from [this](https://doi.org/10.1063/1.4928077) publication.

~~~
./cppm-generator --help

USAGE:
    cppm-generator [OPTIONS] --file <FILE>

OPTIONS:
    -b, --bjerrum-length <BJERRUM_LENGTH>    Bjerrum length (Å) [default: 7]
    -h, --help                               Print help information
    -m, --minus <NUM_MINUS>                  Number of negative (-1e) particles [default: 37]
    -N <NUM_TOTAL>                           Total number of particles [default: 643]
    -o, --file <FILE>                        Output structure (.xyz or .pqr)
    -p, --plus <NUM_PLUS>                    Number of positive (+1e) particles [default: 29]
    -r, --radius <RADIUS>                    Sphere radius (Å) [default: 20]
    -s, --steps <STEPS>                      Number of Monte Carlo iterations [default: 10000]
    -V, --version                            Print version information
~~~

## Details and status:

- [x] Random walk on a sphere using spherical coordinates
- [x] Particle-particle interactions using a Coulomb/softcore potential
- [x] Arbitrary mixing of neutral and charged particles
- [x] Output to `.xyz` and `.pqr` files
- [x] Command line interface
- [x] Dipole moment analysis
- [ ] External electric field to induce patches moments
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

Besides scientific use, this project is a first dive
into the Rust programming language.
