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

use crate::particle::Particle;
use itertools::Itertools;

/// Trait for pair energy between two particles
pub trait PairPotential {
    fn energy(&self, particle_1: &Particle, particle_2: &Particle) -> f64;
}

/// Trait for terms in the Hamiltonian (nonbonded etc.)
pub trait EnergyTerm {
    fn energy(&self, particles: &[Particle], indices: &[usize]) -> f64;
}

/// Coulomb interaction + additional soft-core repulsion
pub struct Coulomb {
    /// Bjerrum length, e^2 / 4 x pi x epsilon_0 x epsilon_r * k_B * T
    pub bjerrum_length: f64,
}

impl Coulomb {
    pub fn new(bjerrum_length: f64) -> Self {
        Coulomb { bjerrum_length }
    }
}

impl PairPotential for Coulomb {
    /// Soft repulsive r^12 + Coulomb potential
    fn energy(&self, particle_1: &Particle, particle_2: &Particle) -> f64 {
        let distance = (particle_1.position - particle_2.position).norm();
        4.0 * f64::powi(4.0 / distance, 12)
            + self.bjerrum_length * particle_1.charge * particle_2.charge / distance
    }
}

/// Nonbonded and pair-wise additive interactions
pub struct Nonbonded<T: PairPotential> {
    pair_potential: T,
}

impl<T: PairPotential> Nonbonded<T> {
    pub fn new(pair_potential: T) -> Self {
        Self {
            pair_potential
        }
    }

    /// Sum all pair interactions in vector of particles (kT)
    #[allow(dead_code)]
    pub fn system_energy(&self, particles: &[Particle]) -> f64 {
        let pair_energy = |v: Vec<&Particle>| self.pair_potential.energy(&v[0], &v[1]);
        particles.iter().combinations(2).map(pair_energy).sum::<f64>()
    }

    /// Sum interaction energy of a single particle with all the rest (kT)
    fn particle_energy(&self, particles: &[Particle], index: usize) -> f64 {
        let mut energy = 0.0;
        for (i, particle) in particles.iter().enumerate() {
            if i != index {
                energy += self.pair_potential.energy(particle, &particles[index]);
            }
        }
        energy
    }

    /// Energy of swapping two particles
    fn swap_move_energy(&self, particles: &[Particle], first: usize, second: usize) -> f64 {
        let mut energy: f64 = self.pair_potential.energy(&particles[first], &particles[second]);
        for (i, particle) in particles.iter().enumerate() {
            if i != first && i != second {
                energy += self.pair_potential.energy(&particle, &particles[first])
                    + self.pair_potential.energy(&particle, &particles[second]);
            }
        }
        energy
    }
}

impl<T: PairPotential> EnergyTerm for Nonbonded<T> {
    /// Energy of a subset of particles given by their indices
    fn energy(&self, particles: &[Particle], indices: &[usize]) -> f64 {
        match indices.len() {
            1 => return self.particle_energy(particles, indices[0]),
            2 => return self.swap_move_energy(particles, indices[0], indices[1]),
            _ => panic!("unknown energy request"),
        }
    }
}
