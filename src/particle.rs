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

#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct Particle {
    pub charge: f64,
    pub phi: f64,   // 0 ≤ φ < 2π (ISO standard)
    pub theta: f64, // 0 ≤ θ ≤ π (ISO standard)
    pub radius: f64, // spherical coordindate radius
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
    pub fn set_angles(&mut self, phi : f64, theta : f64) {
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

        let new_phi = self.phi + dp_phi * f64::sin(2.0 * PI * random::<f64>());
        let new_theta = self.theta + dp_theta * f64::acos(2.0 * random::<f64>() - 1.0);

        //let new_phi = 2.0 * PI * random::<f64>();
        // let new_theta = f64::acos(2.0 * random::<f64>() - 1.0);
        self.set_angles(new_phi, new_theta);
    }
}
