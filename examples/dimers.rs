#![allow(non_snake_case)]
use bebop::define_system;

define_system! {
    Dimers { gene, mRNA, protein, dimer }
    r_transcription : gene             => gene, mRNA    @ 25.
    r_translation   : mRNA             => mRNA, protein @ 1000.
    r_dimerization  : protein, protein => dimer         @ 0.001
    r_decay_mRNA    : mRNA             =>               @ 0.1
    r_decay_protein : protein          =>               @ 1.
}

fn main() {
    let mut num = Vec::new();
    for _ in 0..10000 {
        let mut problem = Dimers::new();
        problem.gene = 1;
        problem.advance_until(1.);
        num.push(problem.dimer);
    }
    println!("{}", num.iter().sum::<isize>() as f64 / num.len() as f64);
}
