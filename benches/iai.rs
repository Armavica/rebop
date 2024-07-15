use iai_callgrind::{library_benchmark, library_benchmark_group, main};
use std::hint::black_box;

use rebop::define_system;
use rebop::gillespie::{Gillespie, Rate};

#[allow(unused_variables)]
fn sir_macro(n: isize) -> isize {
    define_system! {
        r_inf r_heal;
        SIR { S, I, R }
        infection   : S + I => 2 I  @ r_inf
        healing     : I     => R    @ r_heal
    }
    let mut sir = SIR::new();
    sir.seed(0);
    sir.r_inf = 0.1 / 10000.;
    sir.r_heal = 0.05;
    sir.S = n - 1;
    sir.I = 1;
    sir.advance_until(1000.);
    sir.R
}

fn sir_api(n: isize) -> isize {
    let mut sir = Gillespie::new_with_seed([n - 1, 1, 0], 0);
    sir.add_reaction(Rate::lma(0.1 / 10000., [1, 1, 0]), [-1, 1, 0]);
    sir.add_reaction(Rate::lma(0.05, [0, 1, 0]), [0, -1, 1]);
    sir.advance_until(1000.);
    sir.get_species(1)
}

#[library_benchmark]
#[benches::n(1000, 10000, 100000)]
fn bench_sir_macro(n: isize) -> isize {
    black_box(sir_macro(black_box(n)))
}

#[library_benchmark]
#[benches::n(1000, 10000, 100000)]
fn bench_sir_api(n: isize) -> isize {
    black_box(sir_api(black_box(n)))
}

library_benchmark_group!(
    name = sir_group;
    benchmarks = bench_sir_macro, bench_sir_api
);

main!(library_benchmark_groups = sir_group);
