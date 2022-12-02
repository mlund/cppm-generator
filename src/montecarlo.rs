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

#[cfg(test)]
use crate::num_traits::Float;

use average::Estimate;
use itertools::Itertools;
use rand::prelude::IteratorRandom;
use rand::prelude::SliceRandom;
use rand::random;
use rand::rngs::ThreadRng;
use rand::Rng;

use crate::energy::EnergyTerm;
use crate::particle::Particle;

///
/// Use the Metropolis-Hastings criterion to determine if a
/// move should be accepted or rejected based in the energy difference.
///
/// # Arguments
///
/// * `energy_change` - New energy minus old energy in units of kT
///
fn accept_move(energy_change: f64) -> bool {
    let acceptance_probability = f64::min(1.0, f64::exp(-energy_change));
    random::<f64>() < acceptance_probability
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_accept_move() {
        let max_exponent = f64::ln(f64::max_value());
        assert!(accept_move(-1.0));
        assert!(accept_move(0.0));
        assert!(!accept_move(max_exponent));
        assert!(!accept_move(max_exponent * 1.1));
    }
}

///
/// Interface for Monte Carlo move algorithms that all
/// move schemes should implement.
///
pub trait MoveAlgorithm {
    /// Perform a Metropolis-Hastings Monte Carlo move; returns true if the move was successful.
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool;
}

///
/// Final Monte Carlo move that in addition to a move algorithm, also track
/// acceptance statistics. Instances of `MonteCarloMove` is normally created
/// by `Propagator`
///
struct MonteCarloMove {
    acceptance_ratio: average::Mean,
    move_algorithm: Box<dyn MoveAlgorithm>,
}

impl MonteCarloMove {
    pub fn new(move_algorithm: Box<dyn MoveAlgorithm>) -> Self {
        MonteCarloMove {
            acceptance_ratio: average::Mean::new(),
            move_algorithm,
        }
    }
    /// Ratio of accepted vs. total Monte Carlo moves
    pub fn mean_acceptance(&self) -> f64 {
        self.acceptance_ratio.mean()
    }
}

impl MoveAlgorithm for MonteCarloMove {
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool {
        let accepted = self.move_algorithm.do_move(hamiltonian, particles, rng);
        self.acceptance_ratio.add(accepted as usize as f64);
        accepted
    }
}
///
/// Aggregator for multiple Monte Carlo moves
///
#[derive(Default)]
pub struct Propagator {
    moves: Vec<MonteCarloMove>,
}

impl Propagator {
    // see also here: https://stackoverflow.com/questions/71900568/returning-mutable-reference-of-trait-in-vector
    pub fn push<T: 'static + MoveAlgorithm>(&mut self, move_algorithm: T) {
        self.moves
            .push(MonteCarloMove::new(Box::new(move_algorithm)));
    }

    pub fn print(&self) {
        for (i, _move) in self.moves.iter().enumerate() {
            println!(
                "move {} acceptance ratio = {:.2}",
                i,
                _move.mean_acceptance()
            );
        }
    }
}

impl MoveAlgorithm for Propagator {
    ///
    /// Run randomly selected move
    ///
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool {
        let random_move = self.moves.choose_mut(rng).unwrap();
        random_move.do_move(hamiltonian, particles, rng)
    }
}

///
/// Randomly displace spherical coordinates of a single particle
/// The move is done on a unit disc around the old position
///
#[derive(Builder)]
pub struct DisplaceParticle {
    #[builder(default = "0.01")]
    angular_displacement: f64,
}

impl MoveAlgorithm for DisplaceParticle {
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool {
        let index = rng.gen_range(0..particles.len());
        let particle_backup = particles[index].to_owned();
        let old_energy = hamiltonian.energy(particles, &[index]);

        particles[index].displace_angle(self.angular_displacement);
        let new_energy = hamiltonian.energy(particles, &[index]);
        let energy_change = new_energy - old_energy;
        if !accept_move(energy_change) {
            particles[index].clone_from(&particle_backup); // restore
            return false;
        }
        true
    }
}

///
/// Monte Carlo move to swap charges between two randomly selected particles
///
#[derive(Default)]
pub struct SwapCharges;

impl SwapCharges {
    ///
    /// Swap charges of two particles given by their indices.
    /// This can alternatively be done with the following unsafe code:
    /// ~~~
    /// unsafe {
    ///     let a : *mut f64 = &mut particles[first].charge;
    ///     let b : *mut f64 = &mut particles[second].charge;
    ///     std::ptr::swap(a, b);
    /// }
    /// ~~~
    ///
    fn swap_charges(particles: &mut [Particle], first: usize, second: usize) {
        let mut charge = particles[second].charge;
        std::mem::swap(&mut particles[first].charge, &mut charge);
        std::mem::swap(&mut particles[second].charge, &mut charge);
    }

    ///
    /// Pick two, random and non-repeating particle indices
    ///
    fn random_indices(number_of_particles: usize, rng: &mut ThreadRng) -> (usize, usize) {
        assert!(number_of_particles >= 2);
        let (first, second) = (0..number_of_particles)
            .choose_multiple(rng, 2)
            .iter()
            .copied()
            .collect_tuple()
            .unwrap();
        assert!(first != second);
        (first, second)
    }
}

impl MoveAlgorithm for SwapCharges {
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool {
        let (first, second) = Self::random_indices(particles.len(), rng);
        if particles[first].charge != particles[second].charge {
            let old_energy = hamiltonian.energy(particles, &[first, second]);
            Self::swap_charges(particles, first, second);
            let new_energy = hamiltonian.energy(particles, &[first, second]);
            let energy_change = new_energy - old_energy;
            if !accept_move(energy_change) {
                Self::swap_charges(particles, first, second); // restore old charges
                return false;
            }
        }
        true
    }
}
