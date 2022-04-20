extern crate num_traits;

mod energy;
mod montecarlo;
mod particle;
mod input;
mod output;
mod analysis;

use clap::Parser;
use std::error::Error;
use analysis::Moments;
use montecarlo::{DisplaceParticle, SwapCharges, BareMove};
use particle::generate_particles;

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

    // main Monte Carlo loop
    for _ in 0..args.steps {
        propagator.do_move(&hamiltonian, &mut particles, &mut rng);
        moments.sample(&particles);
    }
    println!("number of samples = {}", args.steps);
    propagator.print();
    moments.print();

    output::save_coordinates(&args.file, &particles)?;
    Ok(())
}
