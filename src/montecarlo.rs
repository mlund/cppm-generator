use crate::energy::{EnergyTerm, PairPotential};
use crate::particle::Particle;
use rand::{random};
use rand::rngs::ThreadRng;
use rand::Rng;
use rand::prelude::SliceRandom;
use rand::prelude::IteratorRandom;
use itertools::Itertools;

pub trait MonteCarloMove {
    fn do_move(&self, hamiltonian: &dyn EnergyTerm, particles: &mut Vec<Particle>, rng: &mut ThreadRng) -> bool;
}

/// Metropolis-Hastings criterion for accepting / rejecting move
///
/// # Arguments
///
/// * `energy_change` - New energy minus old energy in units of kT
///
pub fn accept_move(energy_change: f64) -> bool {
    let acceptance_probability = f64::min(1.0, f64::exp(-energy_change));
    random::<f64>() < acceptance_probability
}

/// Aggregator for multiple Monte Carlo moves picked by random
#[allow(dead_code)]
pub struct Propagator {
    pub moves: Vec<Box<dyn MonteCarloMove>>,
}

impl Propagator {
    pub fn new() -> Self {
        let empty_moves = Vec::<Box<dyn MonteCarloMove>>::new();
        Propagator { moves: empty_moves }
    }

    pub fn push<T: 'static + MonteCarloMove>(&mut self, mc_move: T) {
        self.moves.push(Box::new(mc_move));
    }

    pub fn propagate(&self, hamiltonian: &dyn EnergyTerm, mut particles: &mut Vec<Particle>, mut rng: &mut ThreadRng) -> bool {
        let random_move = self.moves.choose(&mut rng).unwrap();
        let accept = random_move.do_move(hamiltonian, particles, &mut rng);
        accept
    }
}

pub struct DisplaceParticle {}

impl DisplaceParticle {
    pub fn default() -> Self {
        Self {}
    }
}

impl MonteCarloMove for DisplaceParticle {
    fn do_move(&self, hamiltonian: &dyn EnergyTerm, particles: &mut Vec<Particle>, rng: &mut ThreadRng) -> bool {
        let index = rng.gen_range(0..particles.len());
        let particle_backup = particles[index].clone();
        let old_energy = hamiltonian.energy(&particles, &vec!(index));
        particles[index].displace_angle(0.01);
        let new_energy = hamiltonian.energy(&particles, &vec!(index));
        let energy_change = new_energy - old_energy;
        if !accept_move(energy_change) {
            particles[index].clone_from(&particle_backup); // restore
            return false;
        }
        true
    }
}

pub struct SwapCharges {}

impl SwapCharges {
    pub fn default() -> Self {
        Self {}
    }
    /// Swap charges of two particles given by their indices
    fn swap_charges(&self, particles: &mut Vec<Particle>, first: usize, second: usize) {
        let mut charges = particles.iter_mut().map(|i| &mut i.charge);
        std::mem::swap(&mut charges.nth(first), &mut charges.nth(second));
    }
}

impl MonteCarloMove for SwapCharges {
    fn do_move(&self, hamiltonian: &dyn EnergyTerm, particles: &mut Vec<Particle>, rng: &mut ThreadRng) -> bool {
        // pick two random indices
        let (first, second) = (0..particles.len()).choose_multiple(rng, 2).iter().copied().collect_tuple().unwrap();
        if particles[first].charge == 0.0 && particles[first].charge == 0.0 {
            return false;
        }
        let old_energy = hamiltonian.energy(&particles, &vec!(first, second));
        self.swap_charges(particles, first, second);
        let new_energy = hamiltonian.energy(&particles, &vec!(first, second));
        let energy_change = new_energy - old_energy;

        if !accept_move(energy_change) { // reject and...
            self.swap_charges(particles, first, second); // ...swap back charges
            return false;
        }
        true
    }
}
