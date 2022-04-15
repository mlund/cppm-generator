use crate::particle::Particle;
use std::vec::Vec;

/// Trait for pair energy between two particles
pub trait PairPotential {
    fn energy(&self, particle_1: &Particle, particle_2: &Particle) -> f64;
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
        f64::powi(4.0 / distance, 12)
            + self.bjerrum_length * particle_1.charge * particle_2.charge / distance
    }
}

/// Sum all pair interactions in vector of particles (kT)
#[allow(dead_code)]
pub fn system_energy<T: PairPotential>(pair_potential: &T, particles: &Vec<Particle>) -> f64 {
    let mut energy = 0.0;
    for (i, particle_1) in particles.iter().enumerate() {
        for (j, particle_2) in particles.iter().enumerate() {
            if i > j {
                energy += pair_potential.energy(particle_1, particle_2);
            }
        }
    }
    energy
}

/// Sum interaction energy of a single particle with all the rest (kT)
#[allow(dead_code)]
pub fn particle_energy<T: PairPotential>(
    pair_potential: &T,
    particles: &Vec<Particle>,
    index: usize,
) -> f64 {
    let mut energy = 0.0;
    for (i, particle) in particles.iter().enumerate() {
        if i != index {
            energy += pair_potential.energy(particle, &particles[index]);
        }
    }
    energy
}
