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

mod energy;
mod montecarlo;
mod particle;
mod input;
mod output;
mod analysis;

use indicatif::{ProgressBar, ProgressStyle};
use clap::Parser;
use std::error::Error;
use analysis::Moments;
use montecarlo::{DisplaceParticle, SwapCharges, BareMove};
use particle::generate_particles;
use crate::analysis::print_global_properties;

fn main() -> Result<(), Box<dyn Error>> {
    let args = input::Args::parse();
    let mut rng = rand::thread_rng();

    let mut particles = generate_particles(args.radius, args.num_total,
                                           args.num_plus, args.num_minus);

    let pair_potential = energy::Coulomb::new(args.bjerrum_length);
    let hamiltonian = energy::Nonbonded::new(pair_potential);
    let mut moments = Moments::new();

    let mut propagator = montecarlo::Propagator::default();
    propagator.push(DisplaceParticle::default());
    propagator.push(SwapCharges::default());

    // customise progress bar
    let bar = ProgressBar::new(args.steps as u64);
    bar.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    // main Monte Carlo loop
    for i in 0..args.steps {
        if i % 100 == 0 { bar.inc(100) };
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
