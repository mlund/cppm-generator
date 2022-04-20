# `cppm-generator`

## Overview

This is a small command-line tool to generate spherical
_Charged Patchy Particle Models_ (CPPM) which mimics (bio-)colloidal
particles, _e.g._ globular proteins.
Here's an [example use](https://doi.org/10.1063/1.4928077) case:

![Alt text](https://aip.scitation.org/action/showOpenGraphArticleImage?doi=10.1063/1.4928077&id=images/medium/1.4928077.figures.f14.gif "a title")

CPPMs are generated by placing neutral, positive, and negative particles
on the surface of a sphere, and minimise the (free) energy using
Metropolis-Hastings Monte Carlo sampling.

Details and status:

- [x] Random walk on a sphere using spherical coordinates
- [x] Particle-particle interactions using a Coulomb/softcore potential
- [x] Arbitrary mixing of neutral and charged particles
- [x] Output to `.xyz` and `.pqr` files
- [x] Command line interface
- [x] Dipole moment analysis
- [ ] External electric field to induce dipole moments
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
./target/release/cppm-generator -h
~~~

## Motivation

Besides scientific use, this project is a first dive
into the Rust programming language.
