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

fn anint(x: f64) -> f64 {
    let i = if x > 0.0 { (x + 0.5) as i32 } else { (x - 0.5) as i32 };
    i as f64
}

fn displace_linear(angle : f64, dp : f64) -> f64 {
    angle + (random::<f64>() - 0.5) * dp
}

fn displace_cos(angle : f64, dp : f64) -> f64 {
    let sign = if random::<f64>() < 0.5 { 1.0 } else { -1.0 };
    let cosdp = f64::cos(dp);
    let mut _cos = f64::cos(angle) + (2.0 * random::<f64>() - 1.0) * cosdp;
    _cos = _cos - anint(_cos / 2.0) * 2.0; // see Allen + Tildesley 1st ed., p 133
    f64::acos(_cos) * sign
}

fn displace_sin(angle : f64, dp : f64) -> f64 {
    let sign = if random::<f64>() < 0.5 { 1.0 } else { -1.0 };
    let sindp = f64::sin(dp);
    let mut _sin = f64::sin(angle) + (2.0 * random::<f64>() - 1.0) * sindp;
    _sin = _sin - anint(_sin / 2.0) * 2.0;
    f64::asin(_sin) * sign
}


#[allow(dead_code)]
#[derive(Clone, Debug)]
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
        self.set_angles(2.0 * PI * random::<f64>(), PI * random::<f64>());
    }

    pub fn displace_angle(&mut self, dp_phi: f64, dp_theta: f64) {
        // see https://mathworld.wolfram.com/SpherePointPicking.html
        // https://doi.org/10.1016/j.amc.2019.124670

        let random_sign = if random::<f64>() < 0.5 { 1.0 } else { -1.0 };

        // cos PI / 1.99
        // sin 0.01
        let new_phi= displace_cos(self.phi, PI / 1.99);
        let new_theta = displace_linear(self.theta, 0.1);

        //let new_phi = 2.0 * PI * random::<f64>();
        // let new_theta = f64::acos(2.0 * random::<f64>() - 1.0);
        self.set_angles(new_phi, new_theta);
    }
}

/// Generate particle vector with charged and neutral particles
pub fn generate_particles(radius : f64, num_total: usize, num_plus: usize, num_minus: usize) -> Vec<Particle> {
    let mut particles: Vec<Particle> =
        vec![Particle::new(radius, 0.0); num_total];

    if num_plus + num_minus > num_total {
        panic!("number of charged ions exceeds total number of particles")
    }
    // cations in the front; anions in the back; then random positions:
    particles.iter_mut().take(num_plus).for_each(|i| i.charge = 1.0);
    particles.iter_mut().rev().take(num_minus).for_each(|i| i.charge = -1.0);
    particles.iter_mut().for_each(|i| i.random_angles());
    particles
}