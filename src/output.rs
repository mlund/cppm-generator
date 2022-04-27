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
use std::fs::File;
use std::io::Write;

///
/// Save particles to a coordinate file (xyz, pqr, ...)
///
pub fn save_coordinates(filename: &str, particles: &[Particle]) -> std::io::Result<()> {
    if filename.ends_with(".xyz") {
        save_xyzfile(filename, particles)?;
    } else if filename.ends_with(".pqr") {
        save_pqrfile(filename, particles)?;
    } else {
        panic!("file suffix must be .xyz or .pqr") // @todo generate error instead
    }
    Ok(())
}

///
/// Save in XYZ molecular file format (atom names and positions)
///
fn save_xyzfile(filename: &str, particles: &[Particle]) -> std::io::Result<()> {
    let mut xyzfile = File::create(filename)?;
    writeln!(xyzfile, "{}\ngenerated by cppm-generator", particles.len())?;
    for particle in particles {
        let atom_name = deduce_atom_name(particle);
        writeln!(
            xyzfile,
            "{} {} {} {}",
            atom_name, &particle.position[0], &particle.position[1], &particle.position[2]
        )?;
    }
    Ok(())
}

///
/// Save in PQR molecular file format (names, positions, charges, radii)
///
fn save_pqrfile(filename: &str, particles: &[Particle]) -> std::io::Result<()> {
    let mut pqrfile = File::create(filename)?;
    writeln!(pqrfile, "{}\ngenerated by cppm-generator", particles.len())?;
    for (index, particle) in particles.iter().enumerate() {
        let atom_name = deduce_atom_name(particle);
        writeln!(
            pqrfile,
            "{:6}{:5} {:^4.4}{:1}{:3.3} {:1}{:4}{:1}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}",
            "ATOM",
            index + 1,
            atom_name,
            "A",
            "CPP",
            "A",
            1,
            "0",
            &particle.position[0],
            &particle.position[1],
            &particle.position[2],
            &particle.charge,
            2.0
        )?;
    }
    Ok(())
}

///
/// Deduces atom name from the particle charge
///
fn deduce_atom_name(particle: &Particle) -> &str {
    if particle.charge > 0.0 {
        return "PP"; // "Plus" Particle
    }
    if particle.charge < 0.0 {
        return "MP"; // "Minus" Particle
    }
    "NP" // "Neutral" Particle
}
