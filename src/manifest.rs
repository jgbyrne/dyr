use serde::Deserialize;
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
pub struct Manifest {
    pub data: MData,
    pub calibrations: HashMap<String, (f64, f64)>,
    pub constraints: Option<HashMap<String, Vec<String>>>,
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
                ncats: 4,
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

