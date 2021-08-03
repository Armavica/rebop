#![allow(non_snake_case)]
use rebop::define_system;

define_system! {
    rtx rtl rdi rdm rdp;
    Dimers { gene, mRNA, protein, dimer }
    r_transcription : gene      => gene + mRNA      @ rtx
    r_translation   : mRNA      => mRNA + protein   @ rtl
    r_dimerization  : 2 protein => dimer            @ rdi
    r_decay_mRNA    : mRNA      =>                  @ rdm
    r_decay_protein : protein   =>                  @ rdp
}

fn main() {
    let mut num = Vec::new();
    for _ in 0..10000 {
        let mut problem = Dimers::new();
        problem.rtx = 25.;
        problem.rtl = 1000.;
        problem.rdi = 0.001;
        problem.rdm = 0.1;
        problem.rdp = 1.;
        problem.gene = 1;
        problem.advance_until(1.);
        num.push(problem.dimer);
    }
    println!("{}", num.iter().sum::<isize>() as f64 / num.len() as f64);
}
