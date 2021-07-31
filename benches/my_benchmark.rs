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

fn api_sir(c: &mut Criterion) {
    index_enum! { enum SIR { S, I, R } }
    c.bench_function("api_sir", |b| {
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
        r_decay_mrna    : M    =>      @ 0.1
        r_decay_prot    : P    =>      @ 1.
    }
    c.bench_function("macro_dimers", |b| {
        b.iter(|| {
            let mut dimers = Dimers::new();
            dimers.G = 1;
            dimers.advance_until(6.);
        })
    });
}

fn api_dimers(c: &mut Criterion) {
    index_enum! { enum Dimers { G, M, P, D } }
    c.bench_function("api_dimers", |b| {
        b.iter(|| {
            let mut dimers = Gillespie::new([1, 0, 0, 0]);
            dimers.add_reaction(Rate::new(25., &[SRate::LMA(Dimers::G)]), [0, 1, 0, 0]);
            dimers.add_reaction(Rate::new(1000., &[SRate::LMA(Dimers::M)]), [0, 0, 1, 0]);
            dimers.add_reaction(Rate::new(0.001, &[SRate::LMA2(Dimers::P)]), [0, 0, -2, 1]);
            dimers.add_reaction(Rate::new(0.1, &[SRate::LMA(Dimers::M)]), [0, -1, 0, 0]);
            dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::P)]), [0, 0, -1, 0]);
            dimers.advance_until(6.);
        })
    });
}

fn macro_dimers2(c: &mut Criterion) {
    define_system! {
        Dimers { A, A_A, AA }
        r_decay_monomer   : A    =>       @ 1.
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

fn api_dimers2(c: &mut Criterion) {
    index_enum! { enum Dimers { A, A_A, AA } }
    c.bench_function("api_dimers2", |b| {
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

fn macro_mm(c: &mut Criterion) {
    define_system! {
        MM { E, S, ES, P }
        r_forward : E, S => ES      @ 0.0017
        r_backward: ES   => E, S    @ 0.5
        r_cat     : ES   => P       @ 0.1
    }
    c.bench_function("macro_mm", |b| {
        b.iter(|| {
            let mut mm = MM::new();
            mm.E = 301;
            mm.S = 120;
            mm.advance_until(100.);
        })
    });
}

fn macro_vilar(c: &mut Criterion) {
    const alphaA: f64 = 50.;
    const alphapA: f64 = 500.;
    const alphaR: f64 = 0.01;
    const alphapR: f64 = 50.;
    const betaA: f64 = 50.;
    const betaR: f64 = 5.;
    const deltaMA: f64 = 10.;
    const deltaMR: f64 = 0.5;
    const deltaA: f64 = 1.;
    const deltaR: f64 = 0.2;
    const gammaA: f64 = 1.;
    const gammaR: f64 = 1.;
    const gammaC: f64 = 2.;
    const thetaA: f64 = 50.;
    const thetaR: f64 = 100.;
    define_system! {
        Vilar { Da, Dr, Dpa, Dpr, Ma, Mr, A, R, C }
        r_activation_a      : Da, A => Dpa      @ gammaA
        r_activation_r      : Dr, A => Dpr      @ gammaR
        r_deactivation_a    : Dpa   => Da, A    @ thetaA
        r_deactivation_r    : Dpr   => Dr, A    @ thetaR
        r_transcription_a   : Da    => Da, Ma   @ alphaA
        r_transcription_r   : Dr    => Dr, Mr   @ alphaR
        r_transcription_p_a : Dpa   => Dpa, Ma  @ alphapA
        r_transcription_p_r : Dpr   => Dpr, Mr  @ alphapR
        r_translation_a     : Ma    => Ma, A    @ betaA
        r_translation_r     : Mr    => Mr, R    @ betaR
        r_complexation      : A, R  => C        @ gammaC
        r_decomplexation    : C     => R        @ deltaA
        r_decay_mRNA_a      : Ma    =>          @ deltaMA
        r_decay_mRNA_r      : Mr    =>          @ deltaMR
        r_decay_prot_a      : A     =>          @ deltaA
        r_decay_prot_r      : R     =>          @ deltaR
    }
    c.bench_function("macro_vilar", |b| {
        b.iter(|| {
            let mut vilar = Vilar::new();
            vilar.Da = 1;
            vilar.Dr = 1;
            vilar.A = 10;
            vilar.R = 10;
            vilar.C = 10;
            vilar.advance_until(200.);
        })
    });
}

criterion_group!(benches,
    api_sir, macro_sir,
    api_dimers, macro_dimers,
    api_dimers2, macro_dimers2,
    macro_mm,
    macro_vilar,
    );
criterion_main!(benches);
