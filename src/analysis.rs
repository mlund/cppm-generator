use crate::particle::Particle;
use nalgebra::Vector3;

/// Calculates the geometric center
pub fn geometric_center(particles: &Vec<Particle>) -> Vector3<f64> {
    particles.iter().map(|i| i.position).sum::<Vector3<f64>>() / particles.len() as f64
}

/// Calculates the center of charge
pub fn charge_center(particles: &Vec<Particle>) -> Vector3<f64> {
    let absolute_charge = particles.iter().map(|i| f64::abs(i.charge)).sum::<f64>();
    particles
        .iter()
        .map(|i| f64::abs(i.charge) * i.position)
        .sum::<Vector3<f64>>()
        / absolute_charge
}

/// Dipole moment with origin at (0,0,0)
pub fn dipole_moment(particles: &Vec<Particle>) -> Vector3<f64> {
    particles.iter().map(|i| i.charge * i.position).sum()
}

/// Analyze mean geometric center; charge center; and dipole moment
pub struct Moments {
    number_of_samples: u32,
    geometric_center: nalgebra::Vector3<f64>,
    charge_center: nalgebra::Vector3<f64>,
    dipole_moment: nalgebra::Vector3<f64>,
}

impl Moments {
    pub fn new() -> Self {
        Self {
            number_of_samples: 0,
            geometric_center: Vector3::new(0.0, 0.0, 0.0),
            charge_center: Vector3::new(0.0, 0.0, 0.0),
            dipole_moment: Vector3::new(0.0, 0.0, 0.0),
        }
    }
    pub fn sample(&mut self, particles: &Vec<Particle>) {
        self.geometric_center += geometric_center(particles);
        self.charge_center += charge_center(particles);
        self.dipole_moment += dipole_moment(particles);
        self.number_of_samples += 1;
    }

    pub fn print(&self) {
        let cog = self.geometric_center.transpose() / self.number_of_samples as f64;
        println!("cog = |⟨∑𝐫ᵢ/N⟩|     = {:.1} Å", cog.norm());
        let coc = self.charge_center.transpose() / self.number_of_samples as f64;
        println!("coc = |⟨∑|qᵢ|𝐫ᵢ⟩/N| = {:.1} eÅ", coc.norm());
        let mu = self.dipole_moment.transpose() / self.number_of_samples as f64;
        println!("𝛍 = |⟨∑qᵢ𝐫ᵢ⟩|       = {:.1} eÅ", mu.norm());
    }
}
