extern crate fitzroy;
mod manifest;

use std::fs;

fn main() {
    let manstr = fs::read_to_string("./narrow.toml").unwrap();
    let manifest: manifest::Manifest = toml::from_str(&manstr).unwrap();
    let config = manifest.build_config().unwrap();
    let mut mcmc = fitzroy::MCMC::new(&config);

    println!("Root Log Likelihood: {}", mcmc.log_likelihood());
}
