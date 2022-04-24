// =-=-=-=-= main.rs =-=-=-=-=
// This file implements the front-end of the Dyr system
// :: The program builds a configuration from a user-provided `manifest` file
// :: It uses this to construct a fitzroy `MCMC` instance
// :: It then performs the mainloop of the iterative markovian process
// :: The body of each iteration is performed by the fitzroy `MCMC` instance
// :: At regular intervals, logs are printed, and the chain state is snapshotted
// :: At the conclusion of the run a maximum-clade-credibility tree is output

extern crate fitzroy;
mod manifest;

use std::fs;
use std::env;
use std::process;

use fitzroy::*;

pub struct DyrConfig {
    pub steps: usize,
    pub model: cfg::Configuration,
}

fn die(msg: &str) -> ! {
    eprintln!("{}", msg);
    process::exit(1);
}

fn main() {
    // Load contents of manifest file into a String
    let mut args = env::args();
    if args.len() != 2 {
        die("dyr path/to/manifest.toml");
    }

    let manpath = args.nth(1).unwrap();
    let manstr = match fs::read_to_string(&manpath) {
        Ok(s) => s,
        Err(_) => { die("Could not read manifest file")  },
    };

    // Parse `Manifest` and build `DyrConfig`
    let manifest: manifest::Manifest = toml::from_str(&manstr).unwrap_or_else(
        |e| { eprintln!("Error: {}", e); die("Could not parse manifest") });

    let config = manifest.build_config().unwrap_or_else(|_| die("Could not build config"));

    // Create `MCMC` instance
    let mut mcmc = fitzroy::MCMC::new(&config.model);

    // Output initial state
    println!("=-=-= Begin Dyr Runlog =-=-=\n");

    println!("Initialising Run: '{}' for {} steps\n", &manpath, config.steps);
    
    println!("::: Initial Params :::\n{:?}\n{:?}", &mcmc.params.traits, &mcmc.params.tree.prior);
    if manifest.mcmc.debug {
        let initial_tree = fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data);
        println!("\n::: Initial Tree :::\n{}", initial_tree);
    }

    println!("Initial Root Log Likelihood: {}", mcmc.last_log_likelihood);

    println!("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");

    // `Summary` stores snapshotted phylogenies 
    let mut summary = fitzroy::Summary::blank();

    // MCMC iterative procedure 
    for i in 0..config.steps {

        // Attempt MCMC step
        // :: If the step function returns `false` the system has probably hit a bug and is unrecoverable 
        if !mcmc.step() {
            println!("Step {}: Catastrophic Error - dumping tree", i);
            println!("{}", fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data));
        }

        // Every `print` steps output current chain likelihood
        if i % manifest.mcmc.print == 0 {
            let log_posterior_likelihood = mcmc.log_posterior_likelihood();
            println!("Step {:<10}: {:<30} | {:<30}", i ,&mcmc.last_log_likelihood, log_posterior_likelihood);

            // If `debug` enabled output information about mcmc moves and prior likelihoods
            if manifest.mcmc.debug {
                println!();
                println!("{}", &mcmc.propose.move_ledger());
                println!();
                println!("{}", &mcmc.config.prior_likelihood_ledger(&mcmc.params));
            }
        }

        // Every `dump` steps output current chain parameterisation
        if (i + 1) % manifest.mcmc.dump == 0 {
            println!("::: Interim Params (Step {}) :::\n{:?}\n{:?}", i, &mcmc.params.traits, &mcmc.params.tree.prior);

            let interim_tree = fitzroy::tree::write_newick(&mcmc.params.tree.tree, &config.model.tree.data);
            println!("\n::: Interim Tree (Step {}) :::\n{}", i, interim_tree);
        }

        // After a burn-in period begin taking snapshots of the tree
        if i > 28_000_000 && i % 1000 == 3 {
            summary.snapshot(&mcmc.params.tree.tree);
        }
    }

    // Output final state
    println!("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");

    println!("Concluded Dyr Run '{}' for {} steps\n", &manpath, config.steps);

    println!("Final Root Log Likelihood: {}\n", mcmc.last_log_likelihood);
    println!("::: Final Params :::\n{:?}\n{:?}", &mcmc.params.traits, &mcmc.params.tree.prior);

    let mcc = fitzroy::tree::write_newick(&summary.mcc_tree(), &config.model.tree.data);
    println!("\n::: MCC Tree ({} snapshots) :::\n{}", &summary.snapshots, mcc);

    println!("\n=-=-= End Dyr Runlog =-=-=");
}
