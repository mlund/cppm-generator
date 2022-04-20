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