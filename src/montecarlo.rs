use crate::energy::PairPotential;
use crate::particle::Particle;
use rand::random;
use rand::rngs::ThreadRng;
use rand::Rng;

/// Metropolis-Hastings criterion for accepting / rejecting move
///
/// # Arguments
///
/// * `energy_change` - New energy minus old energy in units of kT
///
pub fn accept_move(energy_change: f64) -> bool {
    let acceptance_probability = f64::min(1.0, f64::exp(-energy_change));
    random::<f64>() < acceptance_probability
}

pub fn propagate<T: PairPotential>(
    pair_potential: &T,
    mut particles: &mut Vec<Particle>,
    mut rng: &mut ThreadRng,
) -> bool {
    displace_random_particle(pair_potential, &mut particles, &mut rng)
}

fn displace_random_particle<T: PairPotential>(
    pair_potential: &T,
    particles: &mut Vec<Particle>,
    rng: &mut &mut ThreadRng,
) -> bool {
    let index = rng.gen_range(0..particles.len());
    let particle_backup = particles[index].clone();
    let old_energy = crate::energy::particle_energy(pair_potential, &particles, index);
    particles[index].displace_angle(0.05, 0.05);
    let new_energy = crate::energy::particle_energy(pair_potential, &particles, index);
    let energy_change = new_energy - old_energy;
    if !accept_move(energy_change) {
        particles[index].clone_from(&particle_backup); // restore
        return false;
    }
    true
}
