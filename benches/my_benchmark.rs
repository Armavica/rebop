#![allow(unused_variables)]
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rebop::define_system;
use rebop::gillespie::{Gillespie, Rate};

fn bench_sir(c: &mut Criterion) {
    define_system! {
        r_inf r_heal;
        SIR { S, I, R }
        infection   : S + I => 2 I  @ r_inf
        healing     : I     => R    @ r_heal
    }
    let mut group = c.benchmark_group("sir");
    for n in &[10_000, 100_000, 1_000_000] {
        group.bench_with_input(BenchmarkId::new("api", n), n, |b, n| {
            b.iter(|| {
                let mut sir = Gillespie::new_with_seed([n - 1, 1, 0], false, 0);
                sir.add_reaction(Rate::lma(0.1 / 10000., [1, 1, 0]), [-1, 1, 0]);
                sir.add_reaction(Rate::lma(0.05, [0, 1, 0]), [0, -1, 1]);
                sir.advance_until(1000.);
            })
        });
        group.bench_with_input(BenchmarkId::new("macro", n), n, |b, n| {
            b.iter(|| {
                let mut sir = SIR::new();
                sir.seed(0);
                sir.r_inf = 0.1 / 10000.;
                sir.r_heal = 0.05;
                sir.S = n - 1;
                sir.I = 1;
                sir.advance_until(1000.);
            })
        });
    }
    group.finish();
}

fn bench_dimers(c: &mut Criterion) {
    define_system! {
        r_tx r_tl r_dim r_decay_mrna r_decay_prot;
        Dimers { G, M, P, D }
        transcription : G    => G + M   @ r_tx
        translation   : M    => M + P   @ r_tl
        dimerization  : 2 P  => D       @ r_dim
        decay_mrna    : M    =>         @ r_decay_mrna
        decay_prot    : P    =>         @ r_decay_prot
    }
    let mut group = c.benchmark_group("dimers");
    group.bench_function("macro", |b| {
        b.iter(|| {
            let mut dimers = Dimers::new();
            dimers.seed(0);
            dimers.r_tx = 25.;
            dimers.r_tl = 1000.;
            dimers.r_dim = 0.001;
            dimers.r_decay_mrna = 0.1;
            dimers.r_decay_prot = 1.;
            dimers.G = 1;
            dimers.advance_until(6.);
        })
    });
    group.bench_function("api", |b| {
        b.iter(|| {
            let mut dimers = Gillespie::new_with_seed([1, 0, 0, 0], false, 0);
            dimers.add_reaction(Rate::lma(25., [1, 0, 0, 0]), [0, 1, 0, 0]);
            dimers.add_reaction(Rate::lma(1000., [0, 1, 0, 0]), [0, 0, 1, 0]);
            dimers.add_reaction(Rate::lma(0.001, [0, 0, 2, 0]), [0, 0, -2, 1]);
            dimers.add_reaction(Rate::lma(0.1, [0, 1, 0, 0]), [0, -1, 0, 0]);
            dimers.add_reaction(Rate::lma(1., [0, 0, 1, 0]), [0, 0, -1, 0]);
            dimers.advance_until(6.);
        })
    });
    group.finish()
}

fn bench_dimers2(c: &mut Criterion) {
    define_system! {
        r_decay_monomer r_dimerization r_monomerization r_irreversible;
        Dimers { A, A_A, AA }
        decay_monomer   : A    =>       @ r_decay_monomer
        dimerization    : 2 A  => A_A   @ r_dimerization
        monomerization  : A_A  => 2 A   @ r_monomerization
        irreversible    : A_A  => AA    @ r_irreversible
    }
    let mut group = c.benchmark_group("dimers2");
    group.bench_function("macro", |b| {
        b.iter(|| {
            let mut dimers = Dimers::new();
            dimers.seed(0);
            dimers.r_decay_monomer = 1.;
            dimers.r_dimerization = 1. / 500.;
            dimers.r_monomerization = 0.5;
            dimers.r_irreversible = 1. / 25.;
            dimers.A = 100000;
            dimers.advance_until(25.);
        })
    });
    group.bench_function("api", |b| {
        b.iter(|| {
            let mut dimers = Gillespie::new_with_seed([100000, 0, 0], false, 0);
            dimers.add_reaction(Rate::lma(1., [1, 0, 0]), [-1, 0, 0]);
            dimers.add_reaction(Rate::lma(1. / 500., [2, 0, 0]), [-2, 1, 0]);
            dimers.add_reaction(Rate::lma(0.5, [0, 1, 0]), [2, -1, 0]);
            dimers.add_reaction(Rate::lma(1. / 25., [0, 1, 0]), [0, -1, 1]);
            dimers.advance_until(25.);
        })
    });
    group.finish()
}

#[rustfmt::skip]
fn api_erk() {
    let k1 = 0.00166;
    let k2 = 0.0001;
    let k3 = 0.1;
    let k4 = 0.00166;
    let h1 = 0.02;
    let h2 = 0.001;
    let h3 = 0.02;
    let h4 = 0.02;
    let km1 = 0.0001;
    let km3 = 0.1;
    let hm1 = 0.02;
    let hm3 = 0.02;
    let mut erk = Gillespie::new_with_seed([201, 100, 10, 20, 2, 23, 10, 10, 10], false, 0);
    // M    MKK     M_MKK   Mp_MKK  MKP     Mpp     Mp  Mpp_MKP     Mp_MKP
    // 1     2        3        4     5       6      7      8          9
    erk.add_reaction(Rate::lma(k1,  [1, 1, 0, 0, 0, 0, 0, 0, 0]), [-1, -1, 1, 0, 0, 0, 0, 0, 0]);
    erk.add_reaction(Rate::lma(km1, [0, 0, 1, 0, 0, 0, 0, 0, 0]), [1, 1, -1, 0, 0, 0, 0, 0, 0]);
    erk.add_reaction(Rate::lma(k2,  [0, 0, 1, 0, 0, 0, 0, 0, 0]), [0, 1, -1, 0, 0, 0, 1, 0, 0]);

    erk.add_reaction(Rate::lma(k3,  [0, 1, 0, 0, 0, 0, 1, 0, 0]), [0, -1, 0, 1, 0, 0, -1, 0, 0]);
    erk.add_reaction(Rate::lma(km3, [0, 0, 0, 1, 0, 0, 0, 0, 0]), [0, 1, 0, -1, 0, 0, 1, 0, 0]);
    erk.add_reaction(Rate::lma(k4,  [0, 0, 0, 1, 0, 0, 0, 0, 0]), [0, 1, 0, -1, 0, 1, 0, 0, 0]);

    erk.add_reaction(Rate::lma(h1,  [0, 0, 0, 0, 1, 1, 0, 0, 0]), [0, 0, 0, 0, -1, -1, 0, 1, 0]);
    erk.add_reaction(Rate::lma(hm1, [0, 0, 0, 0, 0, 0, 0, 1, 0]), [0, 0, 0, 0, 1, 1, 0, -1, 0]);
    erk.add_reaction(Rate::lma(h2,  [0, 0, 0, 0, 0, 0, 0, 1, 0]), [0, 0, 0, 0, 1, 0, 1, -1, 0]);

    erk.add_reaction(Rate::lma(h3,  [0, 0, 0, 0, 1, 0, 1, 0, 0]), [0, 0, 0, 0, -1, 0, -1, 0, 1]);
    erk.add_reaction(Rate::lma(hm3, [0, 0, 0, 0, 0, 0, 0, 0, 1]), [0, 0, 0, 0, 1, 0, 1, 0, -1]);
    erk.add_reaction(Rate::lma(h4,  [0, 0, 0, 0, 0, 0, 0, 0, 1]), [1, 0, 0, 0, 1, 0, 0, 0, -1]);

    erk.advance_until(1000.0);
}

fn macro_erk() {
    define_system! {
        k1 k2 k3 k4 h1 h2 h3 h4 km1 km3 hm1 hm3;
        Erk { M, Mp, Mpp, MKK, MKP, M_MKK, Mp_MKK, Mp_MKP, Mpp_MKP }

        f1: M + MKK     => M_MKK        @ k1
        r1: M_MKK       => M + MKK      @ km1
        i1: M_MKK       => Mp + MKK     @ k2

        f2: Mp + MKK    => Mp_MKK       @ k3
        r2: Mp_MKK      => Mp + MKK     @ km3
        i2: Mp_MKK      => Mpp + MKK    @ k4

        f3: Mpp + MKP   => Mpp_MKP      @ h1
        r3: Mpp_MKP     => Mpp + MKP    @ hm1
        i3: Mpp_MKP     => Mp + MKP     @ h2

        f4: Mp + MKP    => Mp_MKP       @ h3
        r4: Mp_MKP      => Mp + MKP     @ hm3
        i4: Mp_MKP      => M + MKP      @ h4
    }
    let mut erk = Erk::new();
    erk.k1 = 0.00166;
    erk.k2 = 0.0001;
    erk.k3 = 0.1;
    erk.k4 = 0.00166;
    erk.h1 = 0.02;
    erk.h2 = 0.001;
    erk.h3 = 0.02;
    erk.h4 = 0.02;
    erk.km1 = 0.0001;
    erk.km3 = 0.1;
    erk.hm1 = 0.02;
    erk.hm3 = 0.02;
    erk.M = 201;
    erk.MKK = 100;
    erk.M_MKK = 10;
    erk.Mp_MKK = 20;
    erk.MKP = 2;
    erk.Mpp = 23;
    erk.Mp = 10;
    erk.Mpp_MKP = 10;
    erk.Mp_MKP = 10;
    erk.advance_until(1000.0);
}

fn bench_erk(c: &mut Criterion) {
    let mut group = c.benchmark_group("erk");
    group.bench_function("macro", |b| b.iter(|| macro_erk()));
    group.bench_function("api", |b| b.iter(|| api_erk()));
    group.finish();
}

fn bench_mm(c: &mut Criterion) {
    define_system! {
        r_fwd r_bwd r_cat;
        MM { E, S, ES, P }
        forward : E + S => ES       @ r_fwd
        backward: ES    => E + S    @ r_bwd
        cat     : ES    => P        @ r_cat
    }
    let mut group = c.benchmark_group("mm");
    group.bench_function("macro", |b| {
        b.iter(|| {
            let mut mm = MM::new();
            mm.seed(0);
            mm.r_fwd = 0.0017;
            mm.r_bwd = 0.5;
            mm.r_cat = 0.1;
            mm.E = 301;
            mm.S = 120;
            mm.advance_until(100.);
        })
    });
    group.bench_function("api", |b| {
        b.iter(|| {
            let r_fwd = 0.0017;
            let r_bwd = 0.5;
            let r_cat = 0.1;
            let mut mm = Gillespie::new_with_seed([301, 120, 0, 0], false, 0);
            mm.add_reaction(Rate::lma(r_fwd, [1, 1, 0, 0]), [-1, -1, 1, 0]);
            mm.add_reaction(Rate::lma(r_bwd, [0, 0, 1, 0]), [1, 1, -1, 0]);
            mm.add_reaction(Rate::lma(r_cat, [0, 0, 1, 0]), [0, 0, -1, 1]);
            mm.advance_until(100.);
        })
    });
    group.finish();
}

#[rustfmt::skip]
#[allow(non_snake_case)]
fn api_vilar() {
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
    vilar.advance_until(200.);
}

fn bench_vilar(c: &mut Criterion) {
    define_system! {
        alphaA alphapA alphaR alphapR betaA betaR deltaMA deltaMR deltaA deltaR gammaA gammaR gammaC thetaA thetaR;
        Vilar { Da, Dr, Dpa, Dpr, Ma, Mr, A, R, C }
        activation_a        : Da + A    => Dpa      @ gammaA
        activation_r        : Dr + A    => Dpr      @ gammaR
        deactivation_a      : Dpa       => Da + A   @ thetaA
        deactivation_r      : Dpr       => Dr + A   @ thetaR
        transcription_a     : Da        => Da + Ma  @ alphaA
        transcription_r     : Dr        => Dr + Mr  @ alphaR
        transcription_p_a   : Dpa       => Dpa + Ma @ alphapA
        transcription_p_r   : Dpr       => Dpr + Mr @ alphapR
        translation_a       : Ma        => Ma + A   @ betaA
        translation_r       : Mr        => Mr + R   @ betaR
        complexation        : A + R     => C        @ gammaC
        decomplexation      : C         => R        @ deltaA
        decay_mRNA_a        : Ma        =>          @ deltaMA
        decay_mRNA_r        : Mr        =>          @ deltaMR
        decay_prot_a        : A         =>          @ deltaA
        decay_prot_r        : R         =>          @ deltaR
    }
    define_system! {
        alphaA alphapA alphaR alphapR betaA betaR deltaMA deltaMR deltaA deltaR gammaA gammaR gammaC thetaA thetaR;
        VilarBestOrder { Da, Dr, Dpa, Dpr, Ma, Mr, A, R, C }
        translation_a       : Ma        => Ma + A   @ betaA
        complexation        : A + R     => C        @ gammaC
        decomplexation      : C         => R        @ deltaA
        decay_prot_a        : A         =>          @ deltaA
        decay_mRNA_a        : Ma        =>          @ deltaMA
        transcription_p_a   : Dpa       => Dpa + Ma @ alphapA
        translation_r       : Mr        => Mr + R   @ betaR
        decay_prot_r        : R         =>          @ deltaR
        transcription_a     : Da        => Da + Ma  @ alphaA
        activation_r        : Dr + A    => Dpr      @ gammaR
        deactivation_r      : Dpr       => Dr + A   @ thetaR
        activation_a        : Da + A    => Dpa      @ gammaA
        deactivation_a      : Dpa       => Da + A   @ thetaA
        transcription_p_r   : Dpr       => Dpr + Mr @ alphapR
        decay_mRNA_r        : Mr        =>          @ deltaMR
        transcription_r     : Dr        => Dr + Mr  @ alphaR
    }
    define_system! {
        alphaA alphapA alphaR alphapR betaA betaR deltaMA deltaMR deltaA deltaR gammaA gammaR gammaC thetaA thetaR;
        VilarWorstOrder { Da, Dr, Dpa, Dpr, Ma, Mr, A, R, C }
        transcription_r     : Dr        => Dr + Mr  @ alphaR
        decay_mRNA_r        : Mr        =>          @ deltaMR
        transcription_p_r   : Dpr       => Dpr + Mr @ alphapR
        deactivation_a      : Dpa       => Da + A   @ thetaA
        activation_a        : Da + A    => Dpa      @ gammaA
        deactivation_r      : Dpr       => Dr + A   @ thetaR
        activation_r        : Dr + A    => Dpr      @ gammaR
        transcription_a     : Da        => Da + Ma  @ alphaA
        decay_prot_r        : R         =>          @ deltaR
        translation_r       : Mr        => Mr + R   @ betaR
        transcription_p_a   : Dpa       => Dpa + Ma @ alphapA
        decay_mRNA_a        : Ma        =>          @ deltaMA
        decay_prot_a        : A         =>          @ deltaA
        decomplexation      : C         => R        @ deltaA
        complexation        : A + R     => C        @ gammaC
        translation_a       : Ma        => Ma + A   @ betaA
    }
    let mut group = c.benchmark_group("vilar");
    group.bench_function("macro/normal_order", |b| {
        b.iter(|| {
            let mut vilar = VilarBestOrder::with_parameters(
                50., 500., 0.01, 50., 50., 5., 10., 0.5, 1., 0.2, 1., 1., 2., 50., 100.,
            );
            vilar.seed(0);
            vilar.Da = 1;
            vilar.Dr = 1;
            vilar.advance_until(200.);
        })
    });
    group.bench_function("macro/best_order", |b| {
        b.iter(|| {
            let mut vilar = VilarBestOrder::with_parameters(
                50., 500., 0.01, 50., 50., 5., 10., 0.5, 1., 0.2, 1., 1., 2., 50., 100.,
            );
            vilar.seed(0);
            vilar.Da = 1;
            vilar.Dr = 1;
            vilar.advance_until(200.);
        })
    });
    group.bench_function("macro/worst_order", |b| {
        b.iter(|| {
            let mut vilar = VilarWorstOrder::with_parameters(
                50., 500., 0.01, 50., 50., 5., 10., 0.5, 1., 0.2, 1., 1., 2., 50., 100.,
            );
            vilar.seed(0);
            vilar.Da = 1;
            vilar.Dr = 1;
            vilar.advance_until(200.);
        })
    });
    group.bench_function("api/normal_order", |b| b.iter(|| api_vilar()));
    group.finish();
}

fn macro_ring_10() {
    define_system! {
        k;
        Ring { A0, A1, A2, A3, A4, A5, A6, A7, A8, A9 }
        r0  : A0    => A1   @ k
        r1  : A1    => A2   @ k
        r2  : A2    => A3   @ k
        r3  : A3    => A4   @ k
        r4  : A4    => A5   @ k
        r5  : A5    => A6   @ k
        r6  : A6    => A7   @ k
        r7  : A7    => A8   @ k
        r8  : A8    => A9   @ k
        r9  : A9    => A0  @ k
    }
    let mut ring = Ring::new();
    ring.seed(0);
    ring.k = 1.;
    ring.A0 = 1000;
    ring.advance_until(100.);
}

fn macro_ring_20() {
    define_system! {
        k;
        Ring {
            A0, A1, A2, A3, A4, A5, A6, A7, A8, A9,
            A10, A11, A12, A13, A14, A15, A16, A17, A18, A19
        }
        r0  : A0    => A1   @ k
        r1  : A1    => A2   @ k
        r2  : A2    => A3   @ k
        r3  : A3    => A4   @ k
        r4  : A4    => A5   @ k
        r5  : A5    => A6   @ k
        r6  : A6    => A7   @ k
        r7  : A7    => A8   @ k
        r8  : A8    => A9   @ k
        r9  : A9    => A10  @ k
        r10 : A10   => A11  @ k
        r11 : A11   => A12  @ k
        r12 : A12   => A13  @ k
        r13 : A13   => A14  @ k
        r14 : A14   => A15  @ k
        r15 : A15   => A16  @ k
        r16 : A16   => A17  @ k
        r17 : A17   => A18  @ k
        r18 : A18   => A19  @ k
        r19 : A19   => A0   @ k
    }
    let mut ring = Ring::new();
    ring.seed(0);
    ring.k = 1.;
    ring.A0 = 1000;
    ring.advance_until(100.);
}

fn macro_ring_30() {
    define_system! {
        k;
        Ring {
            A0, A1, A2, A3, A4, A5, A6, A7, A8, A9,
            A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,
            A20, A21, A22, A23, A24, A25, A26, A27, A28, A29
        }
        r0  : A0    => A1   @ k
        r1  : A1    => A2   @ k
        r2  : A2    => A3   @ k
        r3  : A3    => A4   @ k
        r4  : A4    => A5   @ k
        r5  : A5    => A6   @ k
        r6  : A6    => A7   @ k
        r7  : A7    => A8   @ k
        r8  : A8    => A9   @ k
        r9  : A9    => A10  @ k
        r10 : A10   => A11  @ k
        r11 : A11   => A12  @ k
        r12 : A12   => A13  @ k
        r13 : A13   => A14  @ k
        r14 : A14   => A15  @ k
        r15 : A15   => A16  @ k
        r16 : A16   => A17  @ k
        r17 : A17   => A18  @ k
        r18 : A18   => A19  @ k
        r19 : A19   => A20  @ k
        r20 : A20   => A21  @ k
        r21 : A21   => A22  @ k
        r22 : A22   => A23  @ k
        r23 : A23   => A24  @ k
        r24 : A24   => A25  @ k
        r25 : A25   => A26  @ k
        r26 : A26   => A27  @ k
        r27 : A27   => A28  @ k
        r28 : A28   => A29  @ k
        r29 : A29   => A0   @ k
    }
    let mut ring = Ring::new();
    ring.seed(0);
    ring.k = 1.;
    ring.A0 = 1000;
    ring.advance_until(100.);
}

fn macro_ring_40() {
    define_system! {
        k;
        Ring {
            A0, A1, A2, A3, A4, A5, A6, A7, A8, A9,
            A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,
            A20, A21, A22, A23, A24, A25, A26, A27, A28, A29,
            A30, A31, A32, A33, A34, A35, A36, A37, A38, A39
        }
        r0  : A0    => A1   @ k
        r1  : A1    => A2   @ k
        r2  : A2    => A3   @ k
        r3  : A3    => A4   @ k
        r4  : A4    => A5   @ k
        r5  : A5    => A6   @ k
        r6  : A6    => A7   @ k
        r7  : A7    => A8   @ k
        r8  : A8    => A9   @ k
        r9  : A9    => A10  @ k
        r10 : A10   => A11  @ k
        r11 : A11   => A12  @ k
        r12 : A12   => A13  @ k
        r13 : A13   => A14  @ k
        r14 : A14   => A15  @ k
        r15 : A15   => A16  @ k
        r16 : A16   => A17  @ k
        r17 : A17   => A18  @ k
        r18 : A18   => A19  @ k
        r19 : A19   => A20  @ k
        r20 : A20   => A21  @ k
        r21 : A21   => A22  @ k
        r22 : A22   => A23  @ k
        r23 : A23   => A24  @ k
        r24 : A24   => A25  @ k
        r25 : A25   => A26  @ k
        r26 : A26   => A27  @ k
        r27 : A27   => A28  @ k
        r28 : A28   => A29  @ k
        r29 : A29   => A30  @ k
        r30 : A30   => A31  @ k
        r31 : A31   => A32  @ k
        r32 : A32   => A33  @ k
        r33 : A33   => A34  @ k
        r34 : A34   => A35  @ k
        r35 : A35   => A36  @ k
        r36 : A36   => A37  @ k
        r37 : A37   => A38  @ k
        r38 : A38   => A39  @ k
        r39 : A39   => A0   @ k
    }
    let mut ring = Ring::new();
    ring.seed(0);
    ring.k = 1.;
    ring.A0 = 1000;
    ring.advance_until(100.);
}

fn macro_ring_50() {
    define_system! {
        k;
        Ring {
            A0, A1, A2, A3, A4, A5, A6, A7, A8, A9,
            A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,
            A20, A21, A22, A23, A24, A25, A26, A27, A28, A29,
            A30, A31, A32, A33, A34, A35, A36, A37, A38, A39,
            A40, A41, A42, A43, A44, A45, A46, A47, A48, A49
        }
        r0  : A0    => A1   @ k
        r1  : A1    => A2   @ k
        r2  : A2    => A3   @ k
        r3  : A3    => A4   @ k
        r4  : A4    => A5   @ k
        r5  : A5    => A6   @ k
        r6  : A6    => A7   @ k
        r7  : A7    => A8   @ k
        r8  : A8    => A9   @ k
        r9  : A9    => A10  @ k
        r10 : A10   => A11  @ k
        r11 : A11   => A12  @ k
        r12 : A12   => A13  @ k
        r13 : A13   => A14  @ k
        r14 : A14   => A15  @ k
        r15 : A15   => A16  @ k
        r16 : A16   => A17  @ k
        r17 : A17   => A18  @ k
        r18 : A18   => A19  @ k
        r19 : A19   => A20  @ k
        r20 : A20   => A21  @ k
        r21 : A21   => A22  @ k
        r22 : A22   => A23  @ k
        r23 : A23   => A24  @ k
        r24 : A24   => A25  @ k
        r25 : A25   => A26  @ k
        r26 : A26   => A27  @ k
        r27 : A27   => A28  @ k
        r28 : A28   => A29  @ k
        r29 : A29   => A30  @ k
        r30 : A30   => A31  @ k
        r31 : A31   => A32  @ k
        r32 : A32   => A33  @ k
        r33 : A33   => A34  @ k
        r34 : A34   => A35  @ k
        r35 : A35   => A36  @ k
        r36 : A36   => A37  @ k
        r37 : A37   => A38  @ k
        r38 : A38   => A39  @ k
        r39 : A39   => A30  @ k
        r40 : A40   => A41  @ k
        r41 : A41   => A42  @ k
        r42 : A42   => A43  @ k
        r43 : A43   => A44  @ k
        r44 : A44   => A45  @ k
        r45 : A45   => A46  @ k
        r46 : A46   => A47  @ k
        r47 : A47   => A48  @ k
        r48 : A48   => A49  @ k
        r49 : A49   => A0   @ k
    }
    let mut ring = Ring::new();
    ring.seed(0);
    ring.k = 1.0;
    ring.A0 = 1000;
    ring.advance_until(100.);
}

fn api_ring(n: usize, k: f64) -> Gillespie {
    let mut x0 = vec![0; n];
    x0[0] = 1000;
    let mut ring = Gillespie::new_with_seed(x0, false, 0);
    for i in 0..n {
        let mut actions = vec![0; n];
        actions[i] -= 1;
        actions[(i + 1) % n] += 1;
        let mut reactants = vec![0; n];
        reactants[i] += 1;
        ring.add_reaction(Rate::lma(k, reactants), actions);
    }
    ring
}

fn bench_ring(c: &mut Criterion) {
    let mut group = c.benchmark_group("ring");
    for n in &[10, 20, 30, 40, 50] {
        group.bench_with_input(BenchmarkId::new("api", n), n, |b, n| {
            b.iter(|| {
                let mut ring = api_ring(*n, 1.0);
                ring.advance_until(100.);
            })
        });
    }
    group.bench_function(BenchmarkId::new("macro", 10), |b| {
        b.iter(|| macro_ring_10())
    });
    group.bench_function(BenchmarkId::new("macro", 20), |b| {
        b.iter(|| macro_ring_20())
    });
    group.bench_function(BenchmarkId::new("macro", 30), |b| {
        b.iter(|| macro_ring_30())
    });
    group.bench_function(BenchmarkId::new("macro", 40), |b| {
        b.iter(|| macro_ring_40())
    });
    group.bench_function(BenchmarkId::new("macro", 50), |b| {
        b.iter(|| macro_ring_50())
    });
    group.finish();
}

fn macro_flocculation_10(x0: isize) {
    define_system! {
        k;
        Flocculation {
            A1, A2, A3, A4, A5, A6, A7, A8, A9, A10
        }
        r0  : 2 A1      => A2 @ k
        r1  : A1 + A2   => A3 @ k
        r2  : A1 + A3   => A4 @ k
        r3  : A1 + A4   => A5 @ k
        r4  : A1 + A5   => A6 @ k
        r5  : A1 + A6   => A7 @ k
        r6  : A1 + A7   => A8 @ k
        r7  : A1 + A8   => A9 @ k
        r8  : A1 + A9   => A10 @ k
        r9  : 2 A2      => A4 @ k
        r10 : A2 + A3   => A5 @ k
        r11 : A2 + A4   => A6 @ k
        r12 : A2 + A5   => A7 @ k
        r13 : A2 + A6   => A8 @ k
        r14 : A2 + A7   => A9 @ k
        r15 : A2 + A8   => A10 @ k
        r16 : 2 A3      => A6 @ k
        r17 : A3 + A4   => A7 @ k
        r18 : A3 + A5   => A8 @ k
        r19 : A3 + A6   => A9 @ k
        r20 : A3 + A7   => A10 @ k
        r21 : 2 A4      => A8 @ k
        r22 : A4 + A5   => A9 @ k
        r23 : A4 + A6   => A10 @ k
        r24 : 2 A5      => A10 @ k
    }
    let mut f = Flocculation::new();
    f.seed(0);
    f.k = 1.;
    f.A1 = x0;
    f.advance_until(1000.);
}

fn api_flocculation(n: usize, k: f64, n0: isize) -> Gillespie {
    let mut x0 = vec![0; n];
    x0[0] = n0;
    let mut flocculation = Gillespie::new_with_seed(x0, false, 0);
    for i in 1..=n / 2 {
        for j in i..=n - i {
            let mut reactants = vec![0; n];
            let mut jump = vec![0; n];
            reactants[i - 1] += 1;
            reactants[j - 1] += 1;
            jump[i - 1] -= 1;
            jump[j - 1] -= 1;
            jump[i + j - 1] += 1;
            flocculation.add_reaction(Rate::lma(k, reactants), jump);
        }
    }
    flocculation
}

fn bench_flocculation(c: &mut Criterion) {
    let mut group = c.benchmark_group("flocculation");
    for x0 in &[1_000, 100_000] {
        for n in &[10, 20, 30, 40, 50] {
            group.bench_with_input(
                BenchmarkId::new("api", format!("{n} {x0}")),
                &(n, x0),
                |b, &(n, x0)| {
                    b.iter(|| {
                        let mut flocculation = api_flocculation(*n, 1.0, *x0);
                        flocculation.advance_until(1000.);
                    })
                },
            );
        }
        group.bench_function(BenchmarkId::new("macro", format!("10 {x0}")), |b| {
            b.iter(|| macro_flocculation_10(*x0))
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_sir,
    bench_dimers,
    bench_dimers2,
    bench_mm,
    bench_erk,
    bench_vilar,
    bench_flocculation,
    bench_ring,
);

criterion_main!(benches);
