#![allow(non_snake_case)]
use bebop::gillespie::*;
use bebop::index_enum;

index_enum! { enum Edda { U1, U2, U7, U8, U5, U6 } }

fn main() {
    for _ in 0..100 {
        let mut problem = Gillespie::new(&[100, 1, 100, 100, 100, 100]);
        let (n1, n2, n3, n4) = (3., 3., 3., 3.);
        let (k1, k2, k3) = (100., 100., 100.);
        let (K1, K2, K3, K4) = (0.6, 0.46, 0.11, 0.11);
        problem.add_reaction(
            Rate::new(
                k1,
                &[
                    SRate::PosHill(Edda::U5, k2 * K1, n1),
                    SRate::NegHill(Edda::U7, k3 * K2, n2),
                ],
            ),
            &[1, 0, 0, 0, 0, 0],
        );
        problem.add_reaction(
            Rate::new(
                k1,
                &[
                    SRate::PosHill(Edda::U6, k2 * K1, n1),
                    SRate::NegHill(Edda::U8, k3 * K2, n2),
                ],
            ),
            &[0, 1, 0, 0, 0, 0],
        );
        problem.add_reaction(
            Rate::new(k3, &[SRate::NegHill(Edda::U1, k1 * K3, n3)]),
            &[0, 0, 1, 0, 0, 0],
        );
        problem.add_reaction(
            Rate::new(k3, &[SRate::NegHill(Edda::U2, k1 * K3, n3)]),
            &[0, 0, 0, 1, 0, 0],
        );
        problem.add_reaction(
            Rate::new(k2, &[SRate::NegHill(Edda::U1, k1 * K4, n4)]),
            &[0, 0, 0, 0, 1, 0],
        );
        problem.add_reaction(
            Rate::new(k2, &[SRate::NegHill(Edda::U2, k1 * K4, n4)]),
            &[0, 0, 0, 0, 0, 1],
        );
        problem.add_reaction(
            Rate::new(1.0, &[SRate::LMA(Edda::U1)]),
            &[-1, 0, 0, 0, 0, 0],
        );
        problem.add_reaction(
            Rate::new(1.0, &[SRate::LMA(Edda::U2)]),
            &[0, -1, 0, 0, 0, 0],
        );
        problem.add_reaction(
            Rate::new(1.0, &[SRate::LMA(Edda::U7)]),
            &[0, 0, -1, 0, 0, 0],
        );
        problem.add_reaction(
            Rate::new(1.0, &[SRate::LMA(Edda::U8)]),
            &[0, 0, 0, -1, 0, 0],
        );
        let _trace = problem.advance_until(100.);
        // println!("{}", trace.iter().map(|(_, s)| s[0]).sum::<isize>());
        // println!("{}", trace.iter().map(|(_, s)| s[1]).sum::<isize>());
        // println!("{}", trace.iter().map(|(_, s)| s[2]).sum::<isize>());
        // println!("{}", trace.iter().map(|(_, s)| s[3]).sum::<isize>());
        // println!("{}", trace.iter().map(|(_, s)| s[4]).sum::<isize>());
        // println!("{}", trace.iter().map(|(_, s)| s[5]).sum::<isize>());
    }
}
