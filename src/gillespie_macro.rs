#[macro_export]
macro_rules! define_system {
    ( $name:ident { $($species:ident),* }
      $($rname:ident : $($r:ident),+ => $($p:ident),+ @ $rate:expr)+
      ) => {
        use rand::distributions::{Uniform, Distribution};
        use rand::thread_rng;
        use rand_distr::Exp;
        /// Structure representing the problem, with the species and the time.
        #[allow(non_snake_case)]
        #[derive(Clone, Debug)]
        struct $name {
            $($species:isize),*,
            t: f64,
        }
        impl $name {
            /// Constructs an object representing the problem, with the
            /// species and the time.
            fn new() -> Self {
                $name {
                    $($species: 0),*,
                    t: 0.,
                }
            }
            /// Proceed and execute the next reaction to happen.
            ///
            /// Returns whether a reaction actually happened (it could
            /// not happen if all the rates are 0).
            #[allow(unused_assignments)]
            fn step(&mut self) -> bool {
                $(let $rname = $crate::reaction! { self: $($r),+ => $($p),+ @ $rate });*;
                let total_rate = 0. $(+ $rname)*;
                if total_rate <= 0. {
                    return false
                }
                let mut rng = thread_rng();
                let exp = Exp::new(total_rate).unwrap();
                self.t += exp.sample(&mut rng);
                let mut reaction_choice = Uniform::new(0., total_rate).sample(&mut rng);
                $crate::choice!(self reaction_choice; $($rname: $($r),* => $($p),*);*);
                true
            }
            /// Simulate the problem until `t = tmax`.
            ///
            /// Returns the list of all changes of the state.
            fn advance_until(&mut self, tmax: f64) -> Vec<Self> {
                let mut r = vec![self.clone()];
                while self.t < tmax {
                    if !self.step() {
                        self.t = tmax;
                    }
                    r.push(self.clone());
                }
                r
            }
        }
    }
}

#[macro_export]
macro_rules! reaction {
    ($self:ident: nil => $($products:ident),* @ $rate:expr) => {
        $rate
    };
    ($self:ident: $($reagents:ident),* => nil @ $rate:expr) => {
        $rate $(* ($self.$reagents as f64))*
    };
    ($self:ident: $($reagents:ident),* => $($products:ident),+ @ $rate:expr) => {
       $rate $(* ($self.$reagents as f64))*
    };
}

#[macro_export]
macro_rules! choice {
    ($self:ident $rc:ident; $rname:ident: nil => $($p:ident),*;
        $($tail:ident: $($tailr:ident),* => $($tailp:ident),*);*) => {
        if $rc < $rname {
            $($self.$p += 1);*
        } else {
            $rc -= $rname;
            $crate::choice!($self $rc; $($tail: $($tailr),* => $($tailp),*);*);
        }
    };
    ($self:ident $rc:ident; $rname:ident: $($r:ident),* => nil;
        $($tail:ident: $($tailr:ident),* => $($tailp:ident),*);*) => {
        if $rc < $rname {
            $($self.$r -= 1);*;
        } else {
            $rc -= $rname;
            $crate::choice!($self $rc; $($tail: $($tailr),* => $($tailp),*);*);
        }
    };
    ($self:ident $rc:ident; $rname:ident: $($r:ident),* => $($p:ident),*;
        $($tail:ident: $($tailr:ident),* => $($tailp:ident),*);*) => {
        if $rc < $rname {
            $($self.$r -= 1);*;
            $($self.$p += 1);*
        } else {
            $rc -= $rname;
            $crate::choice!($self $rc; $($tail: $($tailr),* => $($tailp),*);*);
        }
    };
    ($self:ident $rc:ident; $rname:ident: nil => $($p:ident),*) => {
        $($self.$p += 1);*
    };
    ($self:ident $rc:ident; $rname:ident: $($r:ident),* => nil) => {
        $($self.$r -= 1);*;
    };
    ($self:ident $rc:ident; $rname:ident: $($r:ident),* => $($p:ident),*) => {
        $($self.$r -= 1);*;
        $($self.$p += 1);*
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn sir() {
        define_system! {
            SIR { S, I, R }
            r_infection: S, I => I, I @ 0.1/10000.
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
            r_dm : mRNA             => nil              @ 0.1
            r_dp : protein          => nil              @ 1.
        }
        let mut dimers = Dimers::new();
        dimers.gene = 1;
        dimers.advance_until(1.);
        assert_eq!(dimers.gene, 1);
        assert!(1000 < dimers.dimer);
        assert!(dimers.dimer < 10000);
    }
}


#[cfg(test)]
mod benchmarks {
    use test::Bencher;
    #[bench]
    fn sir(b: &mut Bencher) {
        define_system! {
            SIR { S, I, R }
            r_infection: S, I => I, I @ 0.1/10000.
            r_remission: I => R @ 0.05
        }
        b.iter(|| {
            let mut sir = SIR::new();
            sir.S = 9999;
            sir.I = 1;
            sir.advance_until(1000.);
        });
    }
    #[bench]
    fn dimers(b: &mut Bencher) {
        define_system! {
            Dimers { gene, mRNA, protein, dimer }
            r_tx : gene             => gene, mRNA       @ 25.
            r_tl : mRNA             => mRNA, protein    @ 1000.
            r_di : protein, protein => dimer            @ 0.001
            r_dm : mRNA             => nil              @ 0.1
            r_dp : protein          => nil              @ 1.
        }
        b.iter(|| {
            let mut dimers = Dimers::new();
            dimers.gene = 1;
            dimers.advance_until(1.);
        });
    }
    #[bench]
    fn dimers2(b: &mut Bencher) {
        define_system! {
            Dimers { A, A_A, AA }
            r_da : A => nil @ 1.
            r_di : A, A => A_A @ 1./500.
            r_ud : A_A => A, A @ 0.5
            r_aa : A_A => AA @ 1./25.
        }
        b.iter(|| {
            let mut dimers = Dimers::new();
            dimers.A = 100000;
            dimers.advance_until(30.);
        });
    }
    #[bench]
    fn schloegl(b: &mut Bencher) {
        define_system! {
            Schloegl { A }
            r_1 : nil => A @ 2200.
            r_2 : A => nil @ 37.5
            r_3 : A, A => A, A, A @ 0.18
            r_4 : A, A, A => A, A @ 0.00025
        }
        b.iter(|| {
            let mut problem = Schloegl::new();
            problem.advance_until(10.);
        });
    }
}
