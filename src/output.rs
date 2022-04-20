use crate::particle::Particle;
use std::fs::File;
use std::io::Write;

/// Save particles to a coordinate file (xyz, pqr, ...)
pub fn save_coordinates(filename : &str, particles: &Vec<Particle>) -> std::io::Result<()> {
    if filename.ends_with(".xyz") {
        save_xyzfile(filename, &particles)?;
    } else if filename.ends_with(".pqr") {
        save_pqrfile(filename, &particles)?;
    } else {
        panic!("file suffix must be .xyz or .pqr") // @todo generate error instead
    }
    Ok(())
}

/// Save in XYZ molecular file format (atom names and positions)
fn save_xyzfile(filename: &str, particles: &Vec<Particle>) -> std::io::Result<()> {
    let mut xyzfile = File::create(filename)?;
    write!(
        xyzfile,
        "{}\ngenerated by cppm-generator\n",
        particles.len()
    )?;
    for particle in particles {
        let mut atom_name = "NP";
        if particle.charge > 0.0 {
            atom_name = "PP";
        } else if particle.charge < 0.0 {
            atom_name = "MP";
        }
        write!(
            xyzfile,
            "{} {} {} {}\n",
            atom_name, &particle.position[0], &particle.position[1], &particle.position[2]
        )?;
    }
    Ok(())
}

/// Save in PQR molecular file format (names, positions, charges, radii)
fn save_pqrfile(filename: &str, particles: &Vec<Particle>) -> std::io::Result<()> {
    let mut pqrfile = File::create(filename)?;
    write!(pqrfile,
           "{}\ngenerated by cppm-generator\n", particles.len())?;
    for (index, particle) in particles.iter().enumerate() {
        let mut atom_name = "NP";
        if particle.charge > 0.0 {
            atom_name = "PP";
        } else if particle.charge < 0.0 {
            atom_name = "MP";
        }
        write!(pqrfile,
               "{:6}{:5} {:^4.4}{:1}{:3.3} {:1}{:4}{:1}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}\n",
               "ATOM", index + 1, atom_name, "A", "CPP", "A", 1, "0",
               &particle.position[0], &particle.position[1], &particle.position[2],
               &particle.charge, 2.0)?;
    }
    Ok(())
}
