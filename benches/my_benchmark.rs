use bebop::define_system;
use bebop::gillespie::{AsIndex, Gillespie, Rate, SRate};
use bebop::index_enum;
use criterion::{criterion_group, criterion_main, Criterion};

fn macro_sir(c: &mut Criterion) {
    define_system! {
        SIR { S, I, R }
        r_inf : S, I => I, I @ 0.1/10000.
        r_heal: I => R @ 0.05
    }
    c.bench_function("macro_sir", |b| {
        b.iter(|| {
            let mut sir = SIR::new();
            sir.S = 9999;
            sir.I = 1;
            sir.advance_until(1000.);
        })
    });
}

fn sir(c: &mut Criterion) {
    index_enum! { enum SIR { S, I, R } }
    c.bench_function("sir", |b| {
        b.iter(|| {
            let mut sir = Gillespie::new([9999, 1, 0]);
            sir.add_reaction(
                Rate::new(0.1 / 10000., &[SRate::LMA(SIR::S), SRate::LMA(SIR::I)]),
                [-1, 1, 0],
            );
            sir.add_reaction(Rate::new(0.05, &[SRate::LMA(SIR::R)]), [0, -1, 1]);
            sir.advance_until(1000.);
        })
    });
}

fn macro_dimers(c: &mut Criterion) {
    define_system! {
        Dimers { G, M, P, D }
        r_transcription : G    => G, M @ 25.
        r_translation   : M    => M, P @ 1000.
        r_dimerization  : P, P => D    @ 0.001
        r_decay_mrna    : M    => nil  @ 0.1
        r_decay_prot    : P    => nil  @ 1.
    }
    c.bench_function("macro_dimers", |b| {
        b.iter(|| {
            let mut dimers = Dimers::new();
            dimers.G = 1;
            dimers.advance_until(7.);
        })
    });
}

fn dimers(c: &mut Criterion) {
    index_enum! { enum Dimers { G, M, P, D } }
    c.bench_function("dimers", |b| {
        b.iter(|| {
            let mut dimers = Gillespie::new([1, 0, 0, 0]);
            dimers.add_reaction(Rate::new(25., &[SRate::LMA(Dimers::G)]), [0, 1, 0, 0]);
            dimers.add_reaction(Rate::new(1000., &[SRate::LMA(Dimers::M)]), [0, 0, 1, 0]);
            dimers.add_reaction(Rate::new(0.001, &[SRate::LMA2(Dimers::P)]), [0, 0, -2, 1]);
            dimers.add_reaction(Rate::new(0.1, &[SRate::LMA(Dimers::M)]), [0, -1, 0, 0]);
            dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::P)]), [0, 0, -1, 0]);
            dimers.advance_until(5.);
        })
    });
}

fn macro_dimers2(c: &mut Criterion) {
    define_system! {
        Dimers { A, A_A, AA }
        r_decay_monomer   : A    => nil   @ 1.
        r_dimerization    : A, A => A_A   @ 1. / 500.
        r_monomerization  : A_A  => A, A  @ 0.5
        r_irreversible    : A_A  => AA    @ 1. / 25.
    }
    c.bench_function("macro_dimers2", |b| {
        b.iter(|| {
            let mut dimers = Dimers::new();
            dimers.A = 100000;
            dimers.advance_until(25.);
        })
    });
}

fn dimers2(c: &mut Criterion) {
    index_enum! { enum Dimers { A, A_A, AA } }
    c.bench_function("dimers2", |b| {
        b.iter(|| {
            let mut dimers = Gillespie::new([100000, 0, 0]);
            dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::A)]), [-1, 0, 0]);
            dimers.add_reaction(Rate::new(1. / 500., &[SRate::LMA2(Dimers::A)]), [-2, 1, 0]);
            dimers.add_reaction(Rate::new(0.5, &[SRate::LMA(Dimers::A_A)]), [2, -1, 0]);
            dimers.add_reaction(Rate::new(1. / 25., &[SRate::LMA(Dimers::A_A)]), [0, -1, 1]);
            dimers.advance_until(25.);
        })
    });
}

criterion_group!(benches, sir, macro_sir, dimers, macro_dimers, dimers2, macro_dimers2);
criterion_main!(benches);
