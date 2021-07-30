#[macro_export]
macro_rules! define_system {
    ( $name:ident { $($species:ident),* }
      $($rname:ident : $($r:ident),* => $($p:ident),* @ $rate:expr)+
      ) => {
        use rand::distributions::Distribution;
        use rand::{Rng, SeedableRng};
        use rand::rngs::SmallRng;
        use rand_distr::Exp;
        /// Structure representing the problem, with the species and the time.
        #[allow(non_snake_case)]
        #[derive(Clone, Debug)]
        struct $name {
            $($species:isize),*,
            t: f64,
            rng: SmallRng,
        }
        impl $name {
            /// Constructs an object representing the problem, with the
            /// species and the time.
            fn new() -> Self {
                $name {
                    $($species: 0),*,
                    t: 0.,
                    rng: SmallRng::from_entropy()
                }
            }
            /// Simulate the problem until `t = tmax`.
            fn advance_until(&mut self, tmax: f64) {
                while self.t < tmax {
                    $(let $rname = $rate $(* (self.$r as f64))*;)*
                    let total_rate = 0. $(+ $rname)*;
                    if !(total_rate > 0.) {
                        self.t = tmax;
                        return
                    }
                    self.t += Exp::new(total_rate).unwrap().sample(&mut self.rng);
                    if self.t > tmax {
                        self.t = tmax;
                        return
                    }
                    let reaction_choice = total_rate * self.rng.gen::<f64>();
                    $crate::choice!(self reaction_choice 0.; $($rname: $($r),* => $($p),*);*);
                }
            }
        }
    }
}

#[macro_export]
macro_rules! choice {
    ($self:ident $rc:ident $carry:expr; $rname:ident: $($r:ident),* => $($p:ident),*;
        $($tail:ident: $($tailr:ident),* => $($tailp:ident),*);*) => {
        if $rc < $carry + $rname {
            $($self.$r -= 1;)*
            $($self.$p += 1;)*
        } else {
            $crate::choice!($self $rc $carry + $rname; $($tail: $($tailr),* => $($tailp),*);*);
        }
    };
    ($self:ident $rc:ident $carry:expr; $rname:ident: $($r:ident),* => $($p:ident),*) => {
        $($self.$r -= 1;)*
        $($self.$p += 1;)*
    };
}

#[cfg(test)]
mod tests {
    #[test]
    fn sir() {
        define_system! {
            SIR { S, I, R }
            r_infection: S, I => I, I @ 0.1 / 10000.
            r_remission: I => R @ 0.01
        }
        let mut sir = SIR::new();
        sir.S = 9999;
        sir.I = 1;
        sir.advance_until(1000.);
        assert_eq!(sir.S + sir.I + sir.R, 10000);
    }
    #[test]
    fn dimers() {
        define_system! {
            Dimers { gene, mRNA, protein, dimer }
            r_tx : gene             => gene, mRNA       @ 25.
            r_tl : mRNA             => mRNA, protein    @ 1000.
            r_di : protein, protein => dimer            @ 0.001
            r_dm : mRNA             =>                  @ 0.1
            r_dp : protein          =>                  @ 1.
        }
        let mut dimers = Dimers::new();
        dimers.gene = 1;
        dimers.advance_until(1.);
        assert_eq!(dimers.gene, 1);
        assert!(1000 < dimers.dimer);
        assert!(dimers.dimer < 10000);
    }
    #[test]
    fn birth_death() {
        define_system! {
            BirthDeath { A }
            birth:      => A    @ 10.
            death:  A   =>      @ 0.1
        }
        let mut birth_death = BirthDeath::new();
        birth_death.advance_until(100.);
        assert!(50 < birth_death.A);
        assert!(birth_death.A < 200);
    }
}

// #[cfg(test)]
// mod benchmarks {
//     use test::Bencher;
//     #[bench]
//     fn schloegl(b: &mut Bencher) {
//         define_system! {
//             Schloegl { A }
//             r_1 :     => A @ 2200.
//             r_2 : A =>     @ 37.5
//             r_3 : A, A => A, A, A @ 0.18
//             r_4 : A, A, A => A, A @ 0.00025
//         }
//         b.iter(|| {
//             let mut problem = Schloegl::new();
//             problem.advance_until(10.);
//         });
//     }
// }
