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

/// Trait for pair energy between two particles
pub trait PairPotential {
    fn energy(&self, particle_1: &Particle, particle_2: &Particle) -> f64;
}

/// Trait for terms in the Hamiltonian (nonbonded etc.)
pub trait EnergyTerm {
    fn energy(&self, particles: &[Particle], indices: &[usize]) -> f64;
}

/// Coulomb interaction + additional soft-core repulsion
#[allow(dead_code)]
pub struct Coulomb {
    /// Bjerrum length, e^2 / 4 x pi x epsilon_0 x epsilon_r * k_B * T
    pub bjerrum_length: f64,
}

#[allow(dead_code)]
impl Coulomb {
    pub fn new(bjerrum_length: f64) -> Self {
        Coulomb { bjerrum_length }
    }
}

impl PairPotential for Coulomb {
    fn energy(&self, particle_1: &Particle, particle_2: &Particle) -> f64 {
        let distance = (particle_1.position - particle_2.position).norm();
        4.0 * f64::powi(4.0 / distance, 12)
            + self.bjerrum_length * particle_1.charge * particle_2.charge / distance
    }
}

/// Nonbonded and pair-wise additive interactions
#[allow(dead_code)]
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
    fn system_energy(&self, particles: &[Particle]) -> f64 {
        let mut energy = 0.0;
        for (i, particle_1) in particles.iter().enumerate() {
            for (j, particle_2) in particles.iter().enumerate() {
                if i > j {
                    energy += self.pair_potential.energy(particle_1, particle_2);
                }
            }
        }
        energy
    }

    /// Sum interaction energy of a single particle with all the rest (kT)
    #[allow(dead_code)]
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
    /// @todo should return error if unknown
    fn energy(&self, particles: &[Particle], indices: &[usize]) -> f64 {
        if indices.len() == 1 {
            return self.particle_energy(particles, indices[0]);
        } else if indices.len() == 2 {
            return self.swap_move_energy(particles, indices[0], indices[1]);
        }
        panic!("unknown energy request encountered");
    }
}
