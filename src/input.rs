use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub struct Args {
    /// Output structure (.xyz or .pqr)
    #[clap(short = 'o', long)]
    pub file: String,

    /// Sphere radius (Å)
    #[clap(short = 'r', long, default_value_t = 20.0)]
    pub radius: f64,

    /// Number of Monte Carlo iterations
    #[clap(short, long, default_value_t = 10000)]
    pub steps: u32,

    /// Total number of particles
    #[clap(short = 'N', default_value_t = 643)]
    pub num_total: usize,

    /// Number of positive (+1e) particles
    #[clap(short = 'p', long = "plus", default_value_t = 29)]
    pub num_plus: usize,

    /// Number of negative (-1e) particles
    #[clap(short = 'm', long = "minus", default_value_t = 37)]
    pub num_minus: usize,

    /// Bjerrum length (Å)
    #[clap(short, long, default_value_t = 7.0)]
    pub bjerrum_length: f64,
}