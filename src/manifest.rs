// =-=-=-=-= manifest.rs =-=-=-=-=
// This file parses manifest files and builds run configurations from them


use serde::Deserialize;
use std::collections::HashMap;
use std::fs;

use fitzroy::*;
use fitzroy::util::PriorDist;

use crate::DyrConfig;

// We use `serde` to automatically deserialise from TOML into a `Manifest` 
// :: We maintain a deliberate distinction between the structures we
// :: deserialise into and `DyrConfig`, our internal config struct

#[derive(Deserialize, Debug)]
pub struct MData {
    pub traits: String,
}

pub type MPrior = (String, Vec<f64>);

#[derive(Deserialize, Debug)]
pub struct MPriors {
    pub tree: MPrior,
}

pub fn is_one_hundred() -> usize { 100 }
pub fn is_ten_thousand() -> usize { 10_000 }
pub fn is_false() -> bool { false }

#[derive(Deserialize, Debug)]
pub struct MChain {
    pub steps: usize,
    #[serde(default = "is_one_hundred")]
    pub print: usize,
    #[serde(default = "is_ten_thousand")]
    pub dump: usize,
    #[serde(default = "is_false")]
    pub debug: bool,
}

#[derive(Deserialize, Debug)]
pub struct MConstraint {
    pub name: String,
    pub clade: Vec<String>,
    pub ancestor: Option<String>,
}

#[derive(Deserialize, Debug)]
pub struct Manifest {
    pub data: MData,
    pub calibrations: HashMap<String, (f64, f64)>,
    #[serde(default)]
    pub constraints: Vec<MConstraint>,
    pub priors: MPriors,
    pub mcmc: MChain,
}

impl Manifest {
    // `Manifest` -> `DyrConfig`
    pub fn build_config(&self) -> Result<DyrConfig, ()> {

        // Load NEXUS data from file specified in manifest 
        let nexstr = match fs::read_to_string(&self.data.traits) {
            Ok(s) => s,
            Err(_) => { eprintln!("Could not read NEXUS file"); return Err(()); },
        };

        // Parse NEXUS data with `fitzroy` parser
        let trait_data = match fitzroy::nexus::parse(&nexstr) {
            Ok(nexus) => { nexus.data.unwrap() },
            Err(msg) => { eprintln!("{}", msg); return Err(()); },
        };

        // Convert binary trait data into tip partial likelihoods
        let tips = match trait_data.binary_into_partials() {
            Some(t) => t,
            None => { eprintln!("Could not build tip partials"); return Err(()); },
        };

        // Build tree prior config; we support Uniform { high, low } and Coalescent { num_intervals }
        let tree_prior = match self.priors.tree.0.as_ref() {
            "Uniform" => {
                let params = &self.priors.tree.1;
                if params.len() != 2 { eprintln!("Uniform prior expects two arguments."); return Err(()) }
                cfg::TreePrior::Uniform { root: PriorDist::Uniform { low: params[0], high: params[1] } }
            },
            "Coalescent" => { 
                let params = &self.priors.tree.1;
                if params.len() != 1 { eprintln!("Coalescent prior expects one argument."); return Err(()) }
                cfg::TreePrior::Coalescent { num_intervals: params[0].round() as usize }
            },
            _ => unimplemented!(),
        };

        // Build mapping of tip names to ids
        let mut name_ids = HashMap::new();
        for (i, (s, _)) in tips.iter().enumerate() {
            name_ids.insert(s, i+1);
        }

        // Build tip calibration config
        let mut calibrations = vec![];
        for (tip, (low, high)) in &self.calibrations {
            match name_ids.get(tip) {
                Some(id) => calibrations.push((*id, cfg::Calibration { low: *low, high: *high })),
                None => { eprintln!("No such tip: {}", tip); return Err(()) },
            }
        }

        // Build clade + ancestry constraint config
        // :: Constraints may include earlier constraints as subsets
        // :: Including an ancestry constraint means including the set of
        // :: both the constraint clade and the ancestral language
        
        let mut built: Vec<cfg::Constraint> = vec![];
        let mut cons_idxs: HashMap<String, usize> = HashMap::new();
        for constraint in self.constraints.iter() {
            let mut tip_ids = vec![];
            
            // Iterate through tags (taxa / constraint names) for this constraint
            for tag in constraint.clade.iter() {
                match name_ids.get(tag) {
                    // This tag is a taxon
                    Some(id) => tip_ids.push(*id),
                    None => { 
                        match cons_idxs.get(tag) {
                            // This tag is an earlier constraint
                            Some(con) => {
                                tip_ids.extend(&built[*con].tips);
                                if let Some(tip) = built[*con].ancestor {
                                    tip_ids.push(tip)
                                }
                            },
                            // This tag is invalid
                            None => {
                                eprintln!("Could not understand: {}", tag);
                                return Err(());
                            }
                        }
                    },
                }
            }

            // Store this constraint's name and index so later constraints can include it
            cons_idxs.insert(constraint.name.clone(), built.len());
            
            // Add this constraint to the configuration
            if let Some(tag) = constraint.ancestor.as_ref() {
                // Ancestors must be taxa, not other constraints
                if let Some(id) = name_ids.get(tag) {
                    built.push(cfg::Constraint::ancestry_constraint(tips.len(), *id, tip_ids));
                }
                else {
                    eprintln!("No such tip: {}", tag);
                    return Err(());
                }
            } 
            else {
                built.push(cfg::Constraint::clade_constraint(tips.len(), tip_ids));
            }
        }

        // The `TreeModel` is the tree prior + calibrations + constraints, and the `TreeData` 
        let tree_model = cfg::TreeModel {
            prior: tree_prior,
            data: tree::TreeData::from_tips(trait_data.num_traits().unwrap(), tips),
            calibrations: calibrations,
            constraints: built,
        };

        // The `TraitsModel` is the priors on substitution model variables
        let traits_model = cfg::TraitsModel {
            // For ascertainment bias correction, unimplemented
            num_traits: PriorDist::Reciprocal, 
             
            subst: cfg::SubstitutionModel::BinaryGTR {
                pi_one: PriorDist::Uniform { low: 0.0, high: 1.0 }
            },
             
            base: PriorDist::Reciprocal,
            //base: PriorDist::Exponential { l: 83_000.0 }, /* rama */
             
            asrv: cfg::ASRV {
                enabled: true,
                shape: PriorDist::Exponential { l: 2.5 },
                // shape: PriorDist::Exponential { l: 1.0 }, /* rama */
                ncats: 4,
            },

            abrv: cfg::ABRV {
                enabled: true,
                shape: PriorDist::Exponential { l: 2.5 },
                // shape: PriorDist::Exponential { l: 200.0 }, /* rama */
            },
        };

        // The built config is the fitzroy `Configuration` and the chain length 
        Ok(DyrConfig {
            model: cfg::Configuration {
                tree: tree_model,
                traits: traits_model,
            },
            steps: self.mcmc.steps,
            print: self.mcmc.print,
            dump : self.mcmc.dump,
            debug: self.mcmc.debug,
        })
    }
}

