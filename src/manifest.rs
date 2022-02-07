use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::fs;

use fitzroy::*;
use fitzroy::util::PriorDist;

use crate::DyrConfig;

#[derive(Deserialize, Debug)]
pub struct MData {
    pub traits: String,
}

pub type MPrior = (String, Vec<f64>);

#[derive(Deserialize, Debug)]
pub struct MPriors {
    pub root_age: MPrior,
}

#[derive(Deserialize, Debug)]
pub struct MChain {
    pub steps: usize,
}

#[derive(Deserialize, Debug)]
pub struct Manifest {
    pub data: MData,
    pub calibrations: HashMap<String, (f64, f64)>,
    pub priors: MPriors,
    pub mcmc: MChain,
}

impl Manifest {
    pub fn build_config(&self) -> Result<DyrConfig, ()> {
        let nexstr = fs::read_to_string(&self.data.traits).unwrap();
        let trait_data = match fitzroy::nexus::parse(&nexstr) {
            Ok(nexus) => { nexus.data.unwrap() },
            Err(msg) => { eprintln!("{}", msg); return Err(()); },
        };

        let tree_prior = match self.priors.root_age.0.as_ref() {
            "Uniform" => {
                let params = &self.priors.root_age.1;
                if params.len() != 2 { eprintln!("Uniform prior expects two arguments."); return Err(()) }
                cfg::TreePrior::Uniform { root: PriorDist::Uniform { low: params[0], high: params[1] } }
            },
            _ => unimplemented!(),
        };

        let tips = trait_data.binary_into_partials().unwrap();

        let mut name_ids = HashMap::new();
        for (i, (s, _)) in tips.iter().enumerate() {
            name_ids.insert(s, i+1);
        }

        let mut calibrations = vec![];
        for (tip, (low, high)) in &self.calibrations {
            match name_ids.get(tip) {
                Some(id) => calibrations.push((*id, cfg::Calibration { low: *low, high: *high })),
                None => { eprintln!("No such tip: {}", tip); return Err(()) },
            }
        }

        let tree_model = cfg::TreeModel {
            prior: tree_prior,
            data: tree::TreeData::from_tips(trait_data.num_traits().unwrap(), tips),
            calibrations: calibrations,
            constraints: vec![],
        };

        let traits_model = cfg::TraitsModel {
            num_traits: PriorDist::Reciprocal,
            subst: cfg::SubstitutionModel::BinaryGTR {
                pi_one: PriorDist::Uniform { low: 0.0, high: 1.0 }
            },
            base: PriorDist::Exponential { l: 10000.0 }, /* rama */
            asrv: cfg::ASRV {
                enabled: true,
                shape: PriorDist::Exponential { l: 2.5 },
            },
            abrv: cfg::ABRV {
                enabled: true,
                shape: PriorDist::Exponential { l: 2.5 },
            },
        };

        Ok(DyrConfig {
            steps: self.mcmc.steps,
            model: cfg::Configuration {
                tree: tree_model,
                traits: traits_model,
            },
        })
    }
}

