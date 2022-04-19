extern crate num_traits;

mod energy;
mod montecarlo;
mod particle;
mod input;
mod output;
mod analysis;

use clap::Parser;
use average;
use average::Estimate;
use std::error::Error;
use analysis::Moments;
use crate::montecarlo::{DisplaceParticle, SwapCharges};
use crate::particle::generate_particles;

fn main() -> Result<(), Box<dyn Error>> {
    let args = input::Args::parse();

    let mut particles = generate_particles(args.radius, args.num_total,
                                           args.num_plus, args.num_minus);

    let pair_potential = energy::Coulomb::new(args.bjerrum_length);
    let hamiltonian = energy::Nonbonded::new(pair_potential);
    let mut rng = rand::thread_rng();
    let mut acceptance_ratio = average::Mean::new();
    let mut moments = Moments::new();

    let mut propagator = montecarlo::Propagator::new();
    propagator.push(DisplaceParticle::default());
    propagator.push(SwapCharges::default());

    // main Monte Carlo loop
    for _ in 0..args.steps {
        let accepted = propagator.propagate(&hamiltonian, &mut particles, &mut rng);
        if accepted {
            acceptance_ratio.add(1.0);
        } else {
            acceptance_ratio.add(0.0)
        }
        moments.sample(&particles);
    }

    println!("number of samples = {}", args.steps);
    println!("fraction accepted = {:.3}", acceptance_ratio.mean());
    moments.print();

    output::save_coordinates(&args.file, &particles)?;
    Ok(())
}
