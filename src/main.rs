extern crate fitzroy;
mod manifest;

use fitzroy::*;

use std::fs;
use std::env;
use std::process;

pub struct DyrConfig {
    pub steps: usize,
    pub model: cfg::Configuration,
}

fn main() {
    let mut args = env::args();
    if args.len() != 2 {
        eprintln!("dyr path/to/config.toml");
        process::exit(1);
    }
    let manstr = match fs::read_to_string(args.nth(1).unwrap()) {
        Ok(s) => s,
        Err(_) => { eprintln!("Could not read config file"); process::exit(1); },
    };

    let manifest: manifest::Manifest = toml::from_str(&manstr).unwrap();
    let config = manifest.build_config().unwrap();
    let mut mcmc = fitzroy::MCMC::new(&config.model);
    println!("{}", fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data));
    let init = mcmc.last_log_likelihood;

    println!("Initial Params: {:?} {:?}", &mcmc.params.traits, &mcmc.params.tree.prior);

    for i in 0..config.steps {
        //println!("Step {:?}:", i);
        if !mcmc.step() {
            println!("{}", fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data));
            println!("Step {}: Instability!", i);
        }

        if i % 100 == 0 {
            let log_posterior_likelihood = mcmc.log_posterior_likelihood();
            println!("Step {}: {} | {}", i ,&mcmc.last_log_likelihood, log_posterior_likelihood);
            mcmc.propose.move_log();
        }
    }
    println!("First Root Log Likelihood: {}", init);
    println!("Last  Root Log Likelihood: {}", mcmc.last_log_likelihood);
    println!("Final Params: {:?} {:?}", &mcmc.params.traits, &mcmc.params.tree.prior);
    println!("{}", fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data));
}
