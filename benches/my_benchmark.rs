use criterion::{criterion_group, criterion_main, Criterion};
use rebop::define_system;
use rebop::gillespie::{AsIndex, Gillespie, Rate, SRate};
use rebop::index_enum;

fn macro_sir_10k(c: &mut Criterion) {
    define_system! {
        r_inf r_heal;
        SIR { S, I, R }
        infection   : S + I => 2 I  @ r_inf
        healing     : I     => R    @ r_heal
    }
    c.bench_function("macro_sir_10k", |b| {
        b.iter(|| {
            let mut sir = SIR::new();
            sir.r_inf = 0.1 / 10000.;
            sir.r_heal = 0.05;
            sir.S = 9999;
            sir.I = 1;
            sir.advance_until(1000.);
        })
    });
}

#[allow(non_snake_case)]
fn macro_sir_1M(c: &mut Criterion) {
    define_system! {
        r_inf r_heal;
        SIR { S, I, R }
        infection   : S + I => 2 I  @ r_inf
        healing     : I     => R    @ r_heal
    }
    c.bench_function("macro_sir_1M", |b| {
        b.iter(|| {
            let mut sir = SIR::new();
            sir.r_inf = 0.1 / 10000.;
            sir.r_heal = 0.05;
            sir.S = 999999;
            sir.I = 1;
            sir.advance_until(1000.);
        })
    });
}

fn api_sir_10k(c: &mut Criterion) {
    index_enum! { enum SIR { S, I, R } }
    c.bench_function("api_sir_10k", |b| {
        b.iter(|| {
            let mut sir = Gillespie::new([9999, 1, 0]);
            sir.add_reaction(
                Rate::new(0.1 / 10000., &[SRate::LMA(SIR::S), SRate::LMA(SIR::I)]),
                [-1, 1, 0],
            );
            sir.add_reaction(Rate::new(0.05, &[SRate::LMA(SIR::I)]), [0, -1, 1]);
            sir.advance_until(1000.);
        })
    });
}

#[allow(non_snake_case)]
fn api_sir_1M(c: &mut Criterion) {
    index_enum! { enum SIR { S, I, R } }
    c.bench_function("api_sir_1M", |b| {
        b.iter(|| {
            let mut sir = Gillespie::new([999999, 1, 0]);
            sir.add_reaction(
                Rate::new(0.1 / 10000., &[SRate::LMA(SIR::S), SRate::LMA(SIR::I)]),
                [-1, 1, 0],
            );
            sir.add_reaction(Rate::new(0.05, &[SRate::LMA(SIR::I)]), [0, -1, 1]);
            sir.advance_until(1000.);
        })
    });
}

fn macro_dimers(c: &mut Criterion) {
    define_system! {
        r_tx r_tl r_dim r_decay_mrna r_decay_prot;
        Dimers { G, M, P, D }
        transcription : G    => G + M   @ r_tx
        translation   : M    => M + P   @ r_tl
        dimerization  : 2 P  => D       @ r_dim
        decay_mrna    : M    =>         @ r_decay_mrna
        decay_prot    : P    =>         @ r_decay_prot
    }
    c.bench_function("macro_dimers", |b| {
        b.iter(|| {
            let mut dimers = Dimers::new();
            dimers.r_tx = 25.;
            dimers.r_tl = 1000.;
            dimers.r_dim = 0.001;
            dimers.r_decay_mrna = 0.1;
            dimers.r_decay_prot = 1.;
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
        r_decay_monomer r_dimerization r_monomerization r_irreversible;
        Dimers { A, A_A, AA }
        decay_monomer   : A    =>       @ r_decay_monomer
        dimerization    : 2 A  => A_A   @ r_dimerization
        monomerization  : A_A  => 2 A   @ r_monomerization
        irreversible    : A_A  => AA    @ r_irreversible
    }
    c.bench_function("macro_dimers2", |b| {
        b.iter(|| {
            let mut dimers = Dimers::new();
            dimers.r_decay_monomer = 1.;
            dimers.r_dimerization = 1. / 500.;
            dimers.r_monomerization = 0.5;
            dimers.r_irreversible = 1. / 25.;
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
        r_fwd r_bwd r_cat;
        MM { E, S, ES, P }
        forward : E + S => ES       @ r_fwd
        backward: ES    => E + S    @ r_bwd
        cat     : ES    => P        @ r_cat
    }
    c.bench_function("macro_mm", |b| {
        b.iter(|| {
            let mut mm = MM::new();
            mm.r_fwd = 0.0017;
            mm.r_bwd = 0.5;
            mm.r_cat = 0.1;
            mm.E = 301;
            mm.S = 120;
            mm.advance_until(100.);
        })
    });
}

#[allow(non_snake_case)]
fn api_vilar(c: &mut Criterion) {
    index_enum! { enum Vilar { Da, Dr, Dpa, Dpr, Ma, Mr, A, R, C } }
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
    c.bench_function("api_vilar", |b| {
        b.iter(|| {
            let mut vilar = Gillespie::new([1, 1, 0, 0, 0, 0, 0, 0, 0]);
            vilar.add_reaction(
                Rate::new(gammaA, &[SRate::LMA(Vilar::Da), SRate::LMA(Vilar::A)]),
                [-1, 0, 1, 0, 0, 0, -1, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(gammaR, &[SRate::LMA(Vilar::Dr), SRate::LMA(Vilar::A)]),
                [0, -1, 0, 1, 0, 0, -1, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(thetaA, &[SRate::LMA(Vilar::Dpa)]),
                [1, 0, -1, 0, 0, 0, 1, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(thetaR, &[SRate::LMA(Vilar::Dpr)]),
                [0, 1, 0, -1, 0, 0, 1, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(alphaA, &[SRate::LMA(Vilar::Da)]),
                [0, 0, 0, 0, 1, 0, 0, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(alphaR, &[SRate::LMA(Vilar::Dr)]),
                [0, 0, 0, 0, 0, 1, 0, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(alphapA, &[SRate::LMA(Vilar::Dpa)]),
                [0, 0, 0, 0, 1, 0, 0, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(alphapR, &[SRate::LMA(Vilar::Dpr)]),
                [0, 0, 0, 0, 0, 1, 0, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(betaA, &[SRate::LMA(Vilar::Ma)]),
                [0, 0, 0, 0, 0, 0, 1, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(betaR, &[SRate::LMA(Vilar::Mr)]),
                [0, 0, 0, 0, 0, 0, 0, 1, 0],
            );
            vilar.add_reaction(
                Rate::new(gammaC, &[SRate::LMA(Vilar::A), SRate::LMA(Vilar::R)]),
                [0, 0, 0, 0, 0, 0, -1, -1, 1],
            );
            vilar.add_reaction(
                Rate::new(gammaA, &[SRate::LMA(Vilar::C)]),
                [0, 0, 0, 0, 0, 0, 0, 1, -1],
            );
            vilar.add_reaction(
                Rate::new(deltaMA, &[SRate::LMA(Vilar::Ma)]),
                [0, 0, 0, 0, -1, 0, 0, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(deltaMR, &[SRate::LMA(Vilar::Mr)]),
                [0, 0, 0, 0, 0, -1, 0, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(deltaA, &[SRate::LMA(Vilar::A)]),
                [0, 0, 0, 0, 0, 0, -1, 0, 0],
            );
            vilar.add_reaction(
                Rate::new(deltaR, &[SRate::LMA(Vilar::R)]),
                [0, 0, 0, 0, 0, 0, 0, -1, 0],
            );
            vilar.advance_until(200.);
        })
    });
}

fn macro_vilar_best_order(c: &mut Criterion) {
    define_system! {
        alphaA alphapA alphaR alphapR betaA betaR deltaMA deltaMR deltaA deltaR gammaA gammaR gammaC thetaA thetaR;
        Vilar { Da, Dr, Dpa, Dpr, Ma, Mr, A, R, C }
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
    c.bench_function("macro_vilar_best_order", |b| {
        b.iter(|| {
            let mut vilar = Vilar::with_parameters(
                50., 500., 0.01, 50., 50., 5., 10., 0.5, 1., 0.2, 1., 1., 2., 50., 100.,
            );
            vilar.Da = 1;
            vilar.Dr = 1;
            vilar.advance_until(200.);
        })
    });
}

fn macro_vilar_worst_order(c: &mut Criterion) {
    define_system! {
        alphaA alphapA alphaR alphapR betaA betaR deltaMA deltaMR deltaA deltaR gammaA gammaR gammaC thetaA thetaR;
        Vilar { Da, Dr, Dpa, Dpr, Ma, Mr, A, R, C }
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
    c.bench_function("macro_vilar_worst_order", |b| {
        b.iter(|| {
            let mut vilar = Vilar::with_parameters(
                50., 500., 0.01, 50., 50., 5., 10., 0.5, 1., 0.2, 1., 1., 2., 50., 100.,
            );
            vilar.Da = 1;
            vilar.Dr = 1;
            vilar.advance_until(200.);
        })
    });
}

fn macro_vilar(c: &mut Criterion) {
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
    c.bench_function("macro_vilar", |b| {
        b.iter(|| {
            let mut vilar = Vilar::with_parameters(
                50., 500., 0.01, 50., 50., 5., 10., 0.5, 1., 0.2, 1., 1., 2., 50., 100.,
            );
            vilar.Da = 1;
            vilar.Dr = 1;
            vilar.advance_until(200.);
        })
    });
}

fn macro_ring_10(c: &mut Criterion) {
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
    c.bench_function("macro_ring_10", |b| {
        b.iter(|| {
            let mut ring = Ring::new();
            ring.k = 1.;
            ring.A0 = 1000;
            ring.advance_until(100.);
        })
    });
}

fn macro_ring_20(c: &mut Criterion) {
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
    c.bench_function("macro_ring_20", |b| {
        b.iter(|| {
            let mut ring = Ring::new();
            ring.k = 1.;
            ring.A0 = 1000;
            ring.advance_until(100.);
        })
    });
}

fn macro_ring_30(c: &mut Criterion) {
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
    c.bench_function("macro_ring_30", |b| {
        b.iter(|| {
            let mut ring = Ring::new();
            ring.k = 1.;
            ring.A0 = 1000;
            ring.advance_until(100.);
        })
    });
}

fn macro_ring_40(c: &mut Criterion) {
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
    c.bench_function("macro_ring_40", |b| {
        b.iter(|| {
            let mut ring = Ring::new();
            ring.k = 1.;
            ring.A0 = 1000;
            ring.advance_until(100.);
        })
    });
}

fn macro_ring_50(c: &mut Criterion) {
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
    c.bench_function("macro_ring_50", |b| {
        b.iter(|| {
            let mut ring = Ring::new();
            ring.k = 1.;
            ring.A0 = 1000;
            ring.advance_until(100.);
        })
    });
}

fn api_ring(n: usize, k: f64) -> Gillespie<usize> {
    let mut x0 = vec![0; n];
    x0[0] = 1000;
    let mut ring = Gillespie::new(x0);
    for i in 0..n {
        let mut actions = vec![0; n];
        actions[i] -= 1;
        actions[(i + 1) % n] += 1;
        ring.add_reaction(Rate::new(k, &[SRate::LMA(i)]), actions);
    }
    ring
}

fn api_ring_10(c: &mut Criterion) {
    c.bench_function("api_ring_10", |b| {
        b.iter(|| {
            let mut ring = api_ring(10, 1.);
            ring.advance_until(100.);
        })
    });
}

fn api_ring_20(c: &mut Criterion) {
    c.bench_function("api_ring_20", |b| {
        b.iter(|| {
            let mut ring = api_ring(20, 1.);
            ring.advance_until(100.);
        })
    });
}

fn api_ring_30(c: &mut Criterion) {
    c.bench_function("api_ring_30", |b| {
        b.iter(|| {
            let mut ring = api_ring(30, 1.);
            ring.advance_until(100.);
        })
    });
}

fn api_ring_40(c: &mut Criterion) {
    c.bench_function("api_ring_40", |b| {
        b.iter(|| {
            let mut ring = api_ring(40, 1.);
            ring.advance_until(100.);
        })
    });
}

fn api_ring_50(c: &mut Criterion) {
    c.bench_function("api_ring_50", |b| {
        b.iter(|| {
            let mut ring = api_ring(50, 1.);
            ring.advance_until(100.);
        })
    });
}

fn macro_flocculation_10(c: &mut Criterion) {
    define_system! {
        k;
        Flocculation {
            A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10
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
    c.bench_function("macro_flocculation_10", |b| {
        b.iter(|| {
            let mut f = Flocculation::new();
            f.k = 1.;
            f.A1 = 1000;
            f.advance_until(100.);
        })
    });
}

criterion_group!(
    benches,
    api_sir_10k,
    macro_sir_10k,
    api_sir_1M,
    macro_sir_1M,
    api_dimers,
    macro_dimers,
    api_dimers2,
    macro_dimers2,
    macro_mm,
    api_vilar,
    macro_vilar,
    macro_vilar_best_order,
    macro_vilar_worst_order,
    macro_ring_10,
    api_ring_10,
    macro_ring_20,
    api_ring_20,
    macro_ring_30,
    api_ring_30,
    macro_ring_40,
    api_ring_40,
    macro_ring_50,
    api_ring_50,
);

criterion_main!(benches);
