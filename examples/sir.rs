#![allow(non_snake_case)]
use rebop::define_system;

define_system! {
    r_inf r_heal;
    SIR { S, I, R }
    infection   : S + I => 2 I @ r_inf
    healing     : I     => R   @ r_heal
}

fn main() {
    let mut num = Vec::new();
    for _ in 0..100000 {
        let mut problem = SIR::new();
        problem.r_inf = 0.1 / 1000.;
        problem.r_heal = 0.01;
        problem.S = 999;
        problem.I = 1;
        problem.advance_until(250.);
        num.push(problem.R);
    }
    println!("{}", num.iter().sum::<isize>() as f64 / num.len() as f64);
}
