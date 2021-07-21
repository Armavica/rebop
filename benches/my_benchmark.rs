use bebop::gillespie::{AsIndex, Gillespie, Rate, SRate};
use bebop::index_enum;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn sir(c: &mut Criterion) {
    index_enum! { enum SIR { S, I, R } }
    c.bench_function("sir", |b| {
        b.iter(|| {
            let mut sir = Gillespie::new(&[9999, 1, 0]);
            sir.add_reaction(
                Rate::new(0.1 / 10000., &[SRate::LMA(SIR::S), SRate::LMA(SIR::I)]),
                &[-1, 1, 0],
            );
            sir.add_reaction(Rate::new(0.05, &[SRate::LMA(SIR::R)]), &[0, -1, 1]);
            sir.advance_until(1000.);
        })
    });
}

fn dimers(c: &mut Criterion) {
    index_enum! { enum Dimers { G, M, P, D } }
    c.bench_function("dimers", |b| {
        b.iter(|| {
            let mut dimers = Gillespie::new(&[1, 0, 0, 0]);
            dimers.add_reaction(Rate::new(25., &[SRate::LMA(Dimers::G)]), &[0, 1, 0, 0]);
            dimers.add_reaction(Rate::new(1000., &[SRate::LMA(Dimers::M)]), &[0, 0, 1, 0]);
            dimers.add_reaction(Rate::new(0.001, &[SRate::LMA2(Dimers::P)]), &[0, 0, -2, 1]);
            dimers.add_reaction(Rate::new(0.1, &[SRate::LMA(Dimers::M)]), &[0, -1, 0, 0]);
            dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::P)]), &[0, 0, -1, 0]);
            dimers.advance_until(1.);
        })
    });
}

fn dimers2(c: &mut Criterion) {
    index_enum! { enum Dimers { A, A_A, AA } }
    c.bench_function("dimers2", |b| {
        b.iter(|| {
            let mut dimers = Gillespie::new(&[100000, 0, 0]);
            dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::A)]), &[-1, 0, 0]);
            dimers.add_reaction(Rate::new(1. / 500., &[SRate::LMA2(Dimers::A)]), &[-2, 1, 0]);
            dimers.add_reaction(Rate::new(0.5, &[SRate::LMA(Dimers::A_A)]), &[2, -1, 0]);
            dimers.add_reaction(Rate::new(1. / 25., &[SRate::LMA(Dimers::A_A)]), &[0, -1, 1]);
            dimers.advance_until(30.);
        })
    });
}

criterion_group!(benches, sir, dimers, dimers2);
criterion_main!(benches);
