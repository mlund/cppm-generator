// Copyright (c) 2022 Mikael Lund
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#[macro_use]
extern crate derive_builder;
extern crate num_traits;

mod analysis;
mod energy;
mod input;
mod montecarlo;
mod output;
mod particle;

use crate::analysis::print_global_properties;
use analysis::Moments;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use montecarlo::{DisplaceParticleBuilder, MonteCarloMove, SwapCharges};
use particle::generate_particles;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let args = input::Args::parse();
    let mut rng = rand::thread_rng();

    // Make particles
    let mut particles =
        generate_particles(args.radius, args.num_total, args.num_plus, args.num_minus);

    // Make Hamiltonian
    let mut hamiltonian = energy::Hamiltonian::default();
    let pair_potential = energy::Coulomb::new(args.bjerrum_length);
    hamiltonian.push(energy::Nonbonded::new(pair_potential));
    if args.target_dipole_moment.is_some() {
        hamiltonian.push(
            energy::ConstrainDipoleBuilder::default()
                .spring_constant(100.0)
                .target_dipole_moment(args.target_dipole_moment.unwrap())
                .build()
                .unwrap(),
        )
    }

    let mut moments = Moments::default();
    let mut propagator = montecarlo::Propagator::default();
    propagator.push(
        DisplaceParticleBuilder::default()
            .angular_displacement(0.01)
            .build()
            .unwrap(),
    );
    propagator.push(SwapCharges::default());

    // customise progress bar
    let bar = ProgressBar::new(args.steps as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
            )
            .progress_chars("#>-"),
    );

    // main Monte Carlo loop
    for i in 0..args.steps {
        if i % 100 == 0 {
            bar.inc(100)
        };
        propagator.do_move(&hamiltonian, &mut particles, &mut rng);
        moments.sample(&particles);
    }
    bar.finish();
    propagator.print();
    moments.print();
    print_global_properties(&particles);

    output::save_coordinates(&args.file, &particles)?;
    Ok(())
}
