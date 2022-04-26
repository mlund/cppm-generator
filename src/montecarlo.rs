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

/// Metropolis-Hastings criterion for accepting / rejecting move
///
/// # Arguments
/// * `energy_change` - New energy minus old energy in units of kT
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
        assert_eq!(accept_move(-1.0), true);
        assert_eq!(accept_move(0.0), true);
        assert_eq!(accept_move(max_exponent), false);
        assert_eq!(accept_move(max_exponent * 1.1), false);
    }
}

/// Trait for Monte Carlo moves
pub trait MonteCarloMove {
    /// Perform a Metropolis-Hastings Monte Carlo move.
    /// Returns true if the move was successful.
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool;
}

/// Properties applicable for all moves
trait MoveProperties {
    /// Average fraction of accepted moves
    fn mean_acceptance(&self) -> f64;
}

/// Trait with both a move function and acceptance statistics
trait MonteCarloMoveExpanded: MonteCarloMove + MoveProperties {}

impl<T: MonteCarloMove + MoveProperties> MonteCarloMoveExpanded for T {}

/// Fully functional MC move with a move function
/// and tracking of acceptance.
struct WrappedMonteCarloMove<T: MonteCarloMove> {
    acceptance_ratio: average::Mean,
    monte_carlo_move: T,
}

impl<T: MonteCarloMove> WrappedMonteCarloMove<T> {
    pub fn new(monte_carlo_move: T) -> Self {
        WrappedMonteCarloMove {
            acceptance_ratio: average::Mean::new(),
            monte_carlo_move,
        }
    }
}

impl<T: MonteCarloMove> MonteCarloMove for WrappedMonteCarloMove<T> {
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool {
        let accepted = self.monte_carlo_move.do_move(hamiltonian, particles, rng);
        self.acceptance_ratio.add(accepted as usize as f64);
        accepted
    }
}

impl<T: MonteCarloMove> MoveProperties for WrappedMonteCarloMove<T> {
    fn mean_acceptance(&self) -> f64 {
        self.acceptance_ratio.mean()
    }
}

/// Aggregator for multiple Monte Carlo moves picked by random
#[derive(Default)]
pub struct Propagator {
    moves: Vec<Box<dyn MonteCarloMoveExpanded>>,
}

impl Propagator {
    // see also here: https://stackoverflow.com/questions/71900568/returning-mutable-reference-of-trait-in-vector
    pub fn push<T: 'static + MonteCarloMove>(&mut self, mc_move: T) {
        let wrapped_move = WrappedMonteCarloMove::new(mc_move);
        self.moves.push(Box::new(wrapped_move));
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

impl MonteCarloMove for Propagator {
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool {
        let random_move = self.moves.choose_mut(rng).unwrap();
        let accepted = random_move.do_move(hamiltonian, particles, rng);
        accepted
    }
}

/// Randomly displace spherical coordinates of a single particle
/// The move is done on a unit disc around the old position
#[derive(Builder)]
pub struct DisplaceParticle {
    #[builder(default = "0.01")]
    angular_displacement: f64,
}

impl MonteCarloMove for DisplaceParticle {
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool {
        let index = rng.gen_range(0..particles.len());
        let particle_backup = particles[index].to_owned();
        let old_energy = hamiltonian.energy(&particles, &[index]);

        particles[index].displace_angle(self.angular_displacement);
        let new_energy = hamiltonian.energy(&particles, &[index]);
        let energy_change = new_energy - old_energy;
        if !accept_move(energy_change) {
            particles[index].clone_from(&particle_backup); // restore
            return false;
        }
        true
    }
}

/// Monte Carlo move to swap charges between two randomly selected particles
#[derive(Default)]
pub struct SwapCharges;

impl SwapCharges {
    /// Swap charges of two particles given by their indices.
    /// Todo list:
    /// - is there a more elegant way to do this using `swap`?
    /// - Mutable charges: `let mut charges = particles.iter_mut().map(|i| &mut i.charge);`
    fn swap_charges(&self, particles: &mut [Particle], first: usize, second: usize) {
        let mut charge = particles[second].charge;
        std::mem::swap(&mut particles[first].charge, &mut charge);
        std::mem::swap(&mut particles[second].charge, &mut charge);
    }
}

impl MonteCarloMove for SwapCharges {
    fn do_move(
        &mut self,
        hamiltonian: &dyn EnergyTerm,
        particles: &mut [Particle],
        rng: &mut ThreadRng,
    ) -> bool {
        // generate two random particle indices
        let (first, second) = (0..particles.len())
            .choose_multiple(rng, 2)
            .iter()
            .copied()
            .collect_tuple()
            .unwrap();
        assert!(first != second);

        if particles[first].charge != particles[second].charge {
            let old_energy = hamiltonian.energy(&particles, &[first, second]);
            self.swap_charges(particles, first, second);
            let new_energy = hamiltonian.energy(&particles, &[first, second]);
            let energy_change = new_energy - old_energy;
            if !accept_move(energy_change) {
                self.swap_charges(particles, first, second); // restore old charges
                return false;
            }
        }
        true
    }
}
