extern crate num_traits;
mod energy;
mod montecarlo;
mod particle;
mod input;
mod output;
mod analysis;

use clap::error::ContextValue::String;
use crate::energy::Coulomb;
use clap::Parser;
use average;
use average::Estimate;
use crate::analysis::Moments;

fn main() {
    let args = input::Args::parse();
    let radius: f64 = args.radius;
    let num_particles = args.num_particles;

    let mut particles: Vec<particle::Particle> =
        vec![particle::Particle::new(radius, 0.0); num_particles];

    if args.num_cations + args.num_anions > args.num_particles {
        // likely wrong way to go...
        panic!("number of charged ions exceeds total number of particles")
    }
    // cations in the front; anions in the back; then random positions.
    particles.iter_mut().take(args.num_cations).for_each(|p| p.charge = 1.0);
    particles.iter_mut().rev().take(args.num_anions).for_each(|p| p.charge = -1.0);
    particles.iter_mut().for_each(|p| p.random_angles());

    let pair_potential = energy::Coulomb::new(args.bjerrum_length);
    let mut rng = rand::thread_rng();
    let mut acceptance_ratio = average::Mean::new();
    let mut moments = Moments::new();

    // main MC loop
    for _ in 0..args.steps {
        let accepted = montecarlo::propagate(&pair_potential, &mut particles, &mut rng);
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

    if args.file.ends_with(".xyz") {
        output::save_xyzfile(&args.file, &particles);
    } else if args.file.ends_with(".pqr") {
        output::save_pqrfile(&args.file, &particles);
    }
}
