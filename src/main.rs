#![allow(non_snake_case)]
#[macro_use]
extern crate bebop;
extern crate rand;

define_system! {
    SIR { S, I, R }
    r_inf : S, I => I, I @ 0.1/10_000.
    r_heal: I => R @ 0.05
}

fn main() {
    for _ in 0..10000 {
        let mut problem = SIR::new();
        problem.S = 9999;
        problem.I = 1;
        let reactions: Vec<SIR> = problem.advance_until(10000.);
        println!("{}", reactions.len());
        // for state in &reactions {
        //     println!("{:?}", state);
        // }
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
