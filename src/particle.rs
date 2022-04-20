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

use nalgebra::Vector3;
use num_traits::Float;
use rand::random;
use std::f64::consts::PI;

/// Convert spherical coordinate to cartesian coordinate
/// See also https://gist.github.com/theypsilon/f09305889e1fd5aa182999af3bad10b9
fn spherical_to_cartesian<T: Float>(phi: T, theta: T, radius: T) -> Vector3<T> {
    Vector3::new(
        radius * phi.sin() * theta.cos(),
        radius * phi.sin() * theta.sin(),
        radius * phi.cos(),
    )
}

/// Particle data incl. position, charge etc.
#[allow(dead_code)]
#[derive(Default, Clone, Debug)]
pub struct Particle {
    pub charge: f64,
    /// 0 ≤ φ < 2π (ISO standard)
    pub phi: f64,
    /// 0 ≤ θ ≤ π (ISO standard)
    pub theta: f64,
    /// spherical coordinate radius
    pub radius: f64,
    /// cartesian position (automatically updated)
    pub position: nalgebra::Vector3<f64>,
}

// associated functions: https://doc.rust-lang.org/rust-by-example/fn/methods.html
#[allow(dead_code)]
impl Particle {
    fn update_cartesian(&mut self) {
        self.position = spherical_to_cartesian(self.phi, self.theta, self.radius);
    }

    pub fn new(radius: f64, charge: f64) -> Self {
        let mut particle = Particle {
            charge: charge,
            phi: 0.0,
            theta: 0.0,
            radius: radius,
            position: Vector3::new(0.0, 0.0, 0.0),
        };
        particle.random_angles();
        particle
    }

    /// Set angles and update cartesian coordinate
    pub fn set_angles(&mut self, phi: f64, theta: f64) {
        self.phi = phi;
        self.theta = theta;
        self.update_cartesian();
    }

    /// Generate random angles and update cartesian coordinate
    pub fn random_angles(&mut self) {
        // see https://mathworld.wolfram.com/SpherePointPicking.html
        let phi = f64::acos(2.0 * random::<f64>() - 1.0);
        let theta = 2.0 * PI * random::<f64>();
        self.set_angles(phi, theta);
    }

    /// Randomly displace theta and phi on a disc
    pub fn displace_angle(&mut self, dp: f64) {
        // see https://mathworld.wolfram.com/SpherePointPicking.html
        // https://doi.org/10.1016/j.amc.2019.124670
        let random_angle = 2.0 * PI * random::<f64>();
        let random_length = dp * random::<f64>();
        let new_phi = self.phi + f64::sin(random_angle) * random_length;
        let new_theta = self.theta + f64::cos(random_angle) * random_length;
        self.set_angles(new_phi, new_theta);
    }
}

/// Generate particle vector with charged and neutral particles
/// @todo Proper error handling / avoid panic
pub fn generate_particles(radius: f64, num_total: usize, num_plus: usize, num_minus: usize) -> Vec<Particle> {
    let mut particles: Vec<Particle> = vec![Particle::new(radius, 0.0); num_total];

    if num_plus + num_minus > num_total {
        panic!("number of charged ions exceeds total number of particles")
    }
    // cations in the front; anions in the back; then random positions:
    particles.iter_mut().take(num_plus).for_each(|i| i.charge = 1.0);
    particles.iter_mut().rev().take(num_minus).for_each(|i| i.charge = -1.0);
    particles.iter_mut().for_each(|i| i.random_angles());
    particles
}