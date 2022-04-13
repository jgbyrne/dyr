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

fn die(msg: &str) -> ! {
    eprintln!("{}", msg);
    process::exit(1);
}

fn main() {
    let mut args = env::args();
    if args.len() != 2 {
        die("dyr path/to/manifest.toml");
    }

    let manpath = args.nth(1).unwrap();
    let manstr = match fs::read_to_string(&manpath) {
        Ok(s) => s,
        Err(_) => { die("Could not read manifest file")  },
    };

    let manifest: manifest::Manifest = toml::from_str(&manstr).unwrap_or_else(
        |e| { eprintln!("Error: {}", e); die("Could not parse manifest") });
    let config = manifest.build_config().unwrap_or_else(|_| die("Could not build config"));
    let mut mcmc = fitzroy::MCMC::new(&config.model);

    println!("=-=-= Begin Dyr Runlog =-=-=\n");

    println!("Initialising Run: '{}' for {} steps\n", &manpath, config.steps);
    
    println!("::: Initial Params :::\n{:?}\n{:?}", &mcmc.params.traits, &mcmc.params.tree.prior);
    if manifest.mcmc.debug {
        let initial_tree = fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data);
        println!("\n::: Initial Tree :::\n{}", initial_tree);
    }

    println!("Initial Root Log Likelihood: {}", mcmc.last_log_likelihood);

    println!("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");

    for i in 0..config.steps {
        //println!("Step {:?}:", i);
        if !mcmc.step() {
            println!("Step {}: Catastrophic Error - dumping tree", i);
            println!("{}", fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data));
        }

        if i % manifest.mcmc.print == 0 {
            let log_posterior_likelihood = mcmc.log_posterior_likelihood();
            println!("Step {:<10}: {:<30} | {:<30}", i ,&mcmc.last_log_likelihood, log_posterior_likelihood);
            if manifest.mcmc.debug {
                mcmc.propose.move_log();
                println!("");
            }
        }

        if (i + 1) % manifest.mcmc.dump == 0 {
            println!("::: Interim Params (Step {}) :::\n{:?}\n{:?}", i, &mcmc.params.traits, &mcmc.params.tree.prior);

            let interim_tree = fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data);
            println!("\n::: Interim Tree (Step {}) :::\n{}", i, interim_tree);
        }

    }

    println!("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");

    println!("Concluded Dyr Run '{}' for {} steps\n", &manpath, config.steps);

    println!("Final Root Log Likelihood: {}\n", mcmc.last_log_likelihood);
    println!("::: Final Params :::\n{:?}\n{:?}", &mcmc.params.traits, &mcmc.params.tree.prior);

    let final_tree = fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data);
    println!("\n::: Final Tree :::\n{}", final_tree);

    println!("\n=-=-= End Dyr Runlog =-=-=");
}
