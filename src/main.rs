extern crate fitzroy;
mod manifest;

use std::fs;

fn main() {
    let manstr = fs::read_to_string("./narrow.toml").unwrap();
    let manifest: manifest::Manifest = toml::from_str(&manstr).unwrap();
    let config = manifest.build_config().unwrap();
    let mut mcmc = fitzroy::MCMC::new(&config);
    println!("{}", fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.tree.data));
    let init = mcmc.last_log_likelihood;

    for i in 0..1000 {
        mcmc.step();
    }
    println!("First Root Log Likelihood: {}", init);
    println!("Last  Root Log Likelihood: {}", mcmc.last_log_likelihood);

}
