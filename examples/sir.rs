#![allow(non_snake_case)]
use bebop::define_system;

define_system! {
    SIR { S, I, R }
    r_inf : S, I => I, I @ 0.1/1000.
    r_heal: I => R @ 0.01
}

fn main() {
    let mut num = Vec::new();
    for _ in 0..100000 {
        let mut problem = SIR::new();
        problem.S = 999;
        problem.I = 1;
        let trace = problem.advance_until(250.);
        num.push(trace.last().unwrap().R);
    }
    println!("{}", num.iter().sum::<isize>() as f64 / num.len() as f64);
}

