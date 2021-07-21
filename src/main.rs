#![allow(non_snake_case)]

use bebop::define_system;

define_system! {
    SIR { S, I, R }
    r_inf : S, I => I, I @ 0.1/1000.
    r_heal: I => R @ 0.01
}

fn main() {
    let mut problem = SIR::new();
    problem.S = 999;
    problem.I = 1;
    let trace: Vec<SIR> = problem.advance_until(250.);
    println!("time,S,I,R");
    for state in &trace {
        println!("{},{},{},{}", state.t, state.S, state.I, state.R);
    }
}

// define_system! {
//     Dimers { gene, mRNA, protein, dimer }
//     r_transcription : gene             => gene, mRNA    @ 25.
//     r_translation   : mRNA             => mRNA, protein @ 1000.
//     r_dimerization  : protein, protein => dimer         @ 0.001
//     r_decay_mRNA    : mRNA             => nil           @ 0.1
//     r_decay_protein : protein          => nil           @ 1.
// }

// fn simulate_once() -> isize {
//     let mut problem = Dimers::new();
//     problem.gene = 1;
//     problem.advance_until(1.);
//     problem.dimer
// }
