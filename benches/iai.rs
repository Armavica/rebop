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
    let mut sir = Gillespie::new_with_seed([n - 1, 1, 0], false, 0);
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

#[rustfmt::skip]
#[allow(non_snake_case)]
fn vilar_api(tmax: f64) {
    // Da, Dr, Dpa, Dpr, Ma, Mr, A, R, C
    let alphaA = 50.;
    let alphapA = 500.;
    let alphaR = 0.01;
    let alphapR = 50.;
    let betaA = 50.;
    let betaR = 5.;
    let deltaMA = 10.;
    let deltaMR = 0.5;
    let deltaA = 1.;
    let deltaR = 0.2;
    let gammaA = 1.;
    let gammaR = 1.;
    let gammaC = 2.;
    let thetaA = 50.;
    let thetaR = 100.;
    let mut vilar = Gillespie::new_with_seed([1, 1, 0, 0, 0, 0, 0, 0, 0], false, 0);
    vilar.add_reaction(Rate::lma(gammaA,  [1, 0, 0, 0, 0, 0, 1, 0, 0]), [-1, 0, 1, 0, 0, 0, -1, 0, 0]);
    vilar.add_reaction(Rate::lma(gammaR,  [0, 1, 0, 0, 0, 0, 1, 0, 0]), [0, -1, 0, 1, 0, 0, -1, 0, 0]);
    vilar.add_reaction(Rate::lma(thetaA,  [0, 0, 1, 0, 0, 0, 0, 0, 0]), [1, 0, -1, 0, 0, 0, 1, 0, 0]);
    vilar.add_reaction(Rate::lma(thetaR,  [0, 0, 0, 1, 0, 0, 0, 0, 0]), [0, 1, 0, -1, 0, 0, 1, 0, 0]);
    vilar.add_reaction(Rate::lma(alphaA,  [1, 0, 0, 0, 0, 0, 0, 0, 0]), [0, 0, 0, 0, 1, 0, 0, 0, 0]);
    vilar.add_reaction(Rate::lma(alphaR,  [0, 1, 0, 0, 0, 0, 0, 0, 0]), [0, 0, 0, 0, 0, 1, 0, 0, 0]);
    vilar.add_reaction(Rate::lma(alphapA, [0, 0, 1, 0, 0, 0, 0, 0, 0]), [0, 0, 0, 0, 1, 0, 0, 0, 0]);
    vilar.add_reaction(Rate::lma(alphapR, [0, 0, 0, 1, 0, 0, 0, 0, 0]), [0, 0, 0, 0, 0, 1, 0, 0, 0]);
    vilar.add_reaction(Rate::lma(betaA,   [0, 0, 0, 0, 1, 0, 0, 0, 0]), [0, 0, 0, 0, 0, 0, 1, 0, 0]);
    vilar.add_reaction(Rate::lma(betaR,   [0, 0, 0, 0, 0, 1, 0, 0, 0]), [0, 0, 0, 0, 0, 0, 0, 1, 0]);
    vilar.add_reaction(Rate::lma(gammaC,  [0, 0, 0, 0, 0, 0, 1, 1, 0]), [0, 0, 0, 0, 0, 0, -1, -1, 1]);
    vilar.add_reaction(Rate::lma(gammaA,  [0, 0, 0, 0, 0, 0, 0, 0, 1]), [0, 0, 0, 0, 0, 0, 0, 1, -1]);
    vilar.add_reaction(Rate::lma(deltaMA, [0, 0, 0, 0, 1, 0, 0, 0, 0]), [0, 0, 0, 0, -1, 0, 0, 0, 0]);
    vilar.add_reaction(Rate::lma(deltaMR, [0, 0, 0, 0, 0, 1, 0, 0, 0]), [0, 0, 0, 0, 0, -1, 0, 0, 0]);
    vilar.add_reaction(Rate::lma(deltaA,  [0, 0, 0, 0, 0, 0, 1, 0, 0]), [0, 0, 0, 0, 0, 0, -1, 0, 0]);
    vilar.add_reaction(Rate::lma(deltaR,  [0, 0, 0, 0, 0, 0, 0, 1, 0]), [0, 0, 0, 0, 0, 0, 0, -1, 0]);
    vilar.advance_until(tmax);
}

#[library_benchmark]
fn bench_vilar_api() {
    black_box(vilar_api(black_box(200.)))
}

library_benchmark_group!(
    name = sir;
    benchmarks = bench_sir_macro, bench_sir_api
);

library_benchmark_group!(
    name = vilar;
    benchmarks = bench_vilar_api
);

main!(library_benchmark_groups = sir, vilar);
