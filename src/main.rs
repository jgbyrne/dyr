extern crate fitzroy;
mod manifest;

use std::fs;

fn main() {
    let manstr = fs::read_to_string("./narrow.toml").unwrap();
    let manifest: manifest::Manifest = toml::from_str(&manstr).unwrap();
    println!("{:#?}", manifest.build_config());

}
