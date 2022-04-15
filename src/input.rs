use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub struct Args {
    #[clap(short = 'o', long)]
    pub file: String,

    /// Sphere radius (Å)
    #[clap(short, long, default_value_t = 20.0)]
    pub radius: f64,

    /// Number of Monte Carlo iterations
    #[clap(short, long, default_value_t = 10000)]
    pub steps: u32,

    /// Total number of particles
    #[clap(short = 'N', long = "nparticles")]
    pub num_particles: usize,

    /// Number of positive (+1e) particles
    #[clap(short = 'p', long = "plus", default_value_t = 20)]
    pub num_cations: usize,

    /// Number of negative (-1e) particles
    #[clap(short = 'm', long = "minus", default_value_t = 29)]
    pub num_anions: usize,

    /// Bjerrum length (Å)
    #[clap(short, long, default_value_t = 7.0)]
    pub bjerrum_length: f64,
}