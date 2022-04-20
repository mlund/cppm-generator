use crate::energy::{EnergyTerm};
use crate::particle::Particle;
use rand::{random};
use rand::rngs::ThreadRng;
use rand::Rng;
use rand::prelude::SliceRandom;
use rand::prelude::IteratorRandom;
use itertools::Itertools;
use average::Estimate;

/// Metropolis-Hastings criterion for accepting / rejecting move
///
/// # Arguments
///
/// * `energy_change` - New energy minus old energy in units of kT
///
fn accept_move(energy_change: f64) -> bool {
    let acceptance_probability = f64::min(1.0, f64::exp(-energy_change));
    random::<f64>() < acceptance_probability
}

/// Trait for Monte Carlo moves
pub trait BareMove {
    fn do_move(&mut self, hamiltonian: &dyn EnergyTerm, particles: &mut Vec<Particle>, rng: &mut ThreadRng) -> bool;
}

/// Properties applicable for all moves
trait MoveProperties {
    fn mean_acceptance(&self) -> f64;
}

/// Trait with both a move function and acceptance statistics
trait MonteCarloMove: BareMove + MoveProperties {}
impl<T: BareMove + MoveProperties> MonteCarloMove for T {}

/// Fully functional MC move with a move function
/// and tracking of acceptance.
struct WrappedMonteCarloMove<T: BareMove> {
    acceptance_ratio: average::Mean,
    monte_carlo_move: T,
}

impl<T: BareMove> WrappedMonteCarloMove<T> {
    pub fn new(monte_carlo_move: T) -> Self {
        WrappedMonteCarloMove {
            acceptance_ratio: average::Mean::new(),
            monte_carlo_move: monte_carlo_move,
        }
    }
}

impl<T: BareMove> BareMove for WrappedMonteCarloMove<T> {
    fn do_move(&mut self, hamiltonian: &dyn EnergyTerm, particles: &mut Vec<Particle>, rng: &mut ThreadRng) -> bool {
        let accepted = self.monte_carlo_move.do_move(hamiltonian, particles, rng);
        self.acceptance_ratio.add(accepted as usize as f64);
        accepted
    }
}

impl<T: BareMove> MoveProperties for WrappedMonteCarloMove<T> {
    fn mean_acceptance(&self) -> f64 {
        self.acceptance_ratio.mean()
    }
}

/// Aggregator for multiple Monte Carlo moves picked by random
#[derive(Default)]
pub struct Propagator {
    moves: Vec<Box<dyn MonteCarloMove>>,
}

impl Propagator {
    // see also here: https://stackoverflow.com/questions/71900568/returning-mutable-reference-of-trait-in-vector
    pub fn push<T: 'static + BareMove>(&mut self, mc_move: T) {
        let wrapped_move = WrappedMonteCarloMove::new(mc_move);
        self.moves.push(Box::new(wrapped_move));
    }

    pub fn print(&self) {
        for (i, _move) in self.moves.iter().enumerate() {
            println!("move {} acceptance ratio = {:.2}", i, _move.mean_acceptance());
        }
    }
}

impl BareMove for Propagator {
    fn do_move(&mut self, hamiltonian: &dyn EnergyTerm, particles: &mut Vec<Particle>, rng: &mut ThreadRng) -> bool {
        let random_move = self.moves.choose_mut(rng).unwrap();
        let accepted = random_move.do_move(hamiltonian, particles, rng);
        accepted
    }
}

/// Randomly displace spherical coordinates of a single particle
/// The move is done on a unit disc around the old position
#[derive(Default)]
pub struct DisplaceParticle;

impl BareMove for DisplaceParticle {
    fn do_move(&mut self, hamiltonian: &dyn EnergyTerm, particles: &mut Vec<Particle>, rng: &mut ThreadRng) -> bool {
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

/// Monte Carlo move to swap charges between two randomly selectec particles
#[derive(Default)]
pub struct SwapCharges;

impl SwapCharges {
    /// Swap charges of two particles given by their indices
    /// @todo is there a more elegant way to do this using `swap`?
    /// Mutable charges: `let mut charges = particles.iter_mut().map(|i| &mut i.charge);`
    fn swap_charges(&self, particles: &mut Vec<Particle>, first: usize, second: usize) {
        let mut charge = particles[second].charge;
        std::mem::swap(&mut particles[first].charge, &mut charge);
        std::mem::swap(&mut particles[second].charge, &mut charge);
    }
}

impl BareMove for SwapCharges {
    fn do_move(&mut self, hamiltonian: &dyn EnergyTerm, particles: &mut Vec<Particle>, rng: &mut ThreadRng) -> bool {
        // pick two random indices
        let (first, second) = (0..particles.len()).choose_multiple(rng, 2).iter().copied().collect_tuple().unwrap();
        assert!(first != second);
        if particles[first].charge != particles[second].charge {
            let old_energy = hamiltonian.energy(&particles, &vec!(first, second));
            self.swap_charges(particles, first, second);
            let new_energy = hamiltonian.energy(&particles, &vec!(first, second));
            let energy_change = new_energy - old_energy;

            if !accept_move(energy_change) { // reject and...
                self.swap_charges(particles, first, second); // ...swap back charges
                return false;
            }
        }
        true
    }
}
