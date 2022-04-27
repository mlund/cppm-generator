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
use nalgebra::Vector3;
use std::f64::consts::PI;

///
/// Total charge
///
fn net_charge(particles: &[Particle]) -> f64 {
    particles.iter().map(|i| i.charge).sum::<f64>()
}

///
/// Total absolute charge
///
fn absolute_charge(particles: &[Particle]) -> f64 {
    particles.iter().map(|i| f64::abs(i.charge)).sum::<f64>()
}

///
/// Calculates the geometric center
///
fn geometric_center(particles: &[Particle]) -> Option<Vector3<f64>> {
    if particles.is_empty() {
        return None;
    }
    Some(particles.iter().map(|i| i.position).sum::<Vector3<f64>>() / particles.len() as f64)
}

///
/// Calculates the center of charge
///
fn charge_center(particles: &[Particle]) -> Vector3<f64> {
    let absolute_charge = particles.iter().map(|i| f64::abs(i.charge)).sum::<f64>();

    particles
        .iter()
        .map(|i| f64::abs(i.charge) * i.position)
        .sum::<Vector3<f64>>()
        / absolute_charge
}

///
/// Dipole moment with origin at (0,0,0)
///
pub fn dipole_moment(particles: &[Particle]) -> Vector3<f64> {
    particles.iter().map(|i| i.charge * i.position).sum()
}

///
/// Analyze mean geometric center; charge center; and dipole moment
///
#[derive(Default)]
pub struct Moments {
    number_of_samples: u32,
    geometric_center: nalgebra::Vector3<f64>,
    charge_center: nalgebra::Vector3<f64>,
    dipole_moment: nalgebra::Vector3<f64>,
    dipole_moment_scalar: f64,
}

impl Moments {
    pub fn sample(&mut self, particles: &[Particle]) {
        self.geometric_center += geometric_center(particles).expect("no particles to sample");
        self.charge_center += charge_center(particles);
        let mu = dipole_moment(particles);
        self.dipole_moment += mu;
        self.dipole_moment_scalar += mu.norm();
        self.number_of_samples += 1;
    }

    pub fn print(&self) {
        let cog = self.geometric_center.transpose() / self.number_of_samples as f64;
        println!(
            "geometric center displacement = |⟨∑𝐫ᵢ/N⟩| = {:.1} Å",
            cog.norm()
        );

        let coc = self.charge_center.transpose() / self.number_of_samples as f64;
        println!(
            "charge center displacement    = |⟨∑|qᵢ|𝐫ᵢ⟩/N| = {:.1} eÅ",
            coc.norm()
        );

        let mu = self.dipole_moment_scalar / self.number_of_samples as f64;
        println!(
            "mean dipole moment 𝛍          = ⟨|∑qᵢ𝐫ᵢ|⟩ = {:.1} eÅ = {:.1} D",
            mu,
            mu / 0.2081943
        );
    }
}

///
/// Print cppm particles such as surface charge density, net charge etc.
///
pub fn print_global_properties(particles: &[Particle]) {
    let radius = particles.first().unwrap().radius;
    let surface_area = 4.0 * PI * radius * radius;
    let mu = dipole_moment(particles).norm();
    println!("CPPM properties:");
    println!("  number of particles       = {}", particles.len());
    println!(
        "  abs. net charge           = {}",
        absolute_charge(particles)
    );
    println!(
        "  radius                    = {} Å",
        particles.first().unwrap().radius
    );
    println!("  surface area              = {:.2} Å²", surface_area);
    println!(
        "  monopole moment           = {:.2}e",
        net_charge(particles)
    );
    println!(
        "  dipole moment |𝛍|         = {:.2} eÅ = {:.2} D",
        mu,
        mu / 0.2081943
    );
    println!(
        "  particle density          = {:.2} Å²/particle",
        surface_area / (particles.len() as f64)
    );
    println!(
        "  surf. charge density      = {:.2} Å²/e",
        surface_area / net_charge(particles)
    );
    println!(
        "  abs. surf. charge density = {:.2} Å²/e",
        surface_area / absolute_charge(particles)
    );
}
