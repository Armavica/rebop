//! Macro-based API to describe chemical reaction networks and simulate
//! them.
//!
//! See [`define_system`].

/// Definition of a chemical reaction network.
///
/// This macro creates a `struct` containing state variables, parameter
/// values and a pseudo-random number generator.  The state variables
/// and parameter values can be modified directly.  It implements three
/// functions: `new`, `with_parameters` and `advance_until`.
///
/// The function `new` creates a new instance of the structure with
/// all state variables set to `0` and all parameter values set to
/// `f64:NAN`.  The parameter values have then to be initialized
/// manually.  If a `NAN` remains at the time of the simulation, no
/// reaction will happen.
///
/// The function `with_parameters` is an alternate initializer that
/// allows to give directly all the parameter values.
///
/// The function `advance_until` simulates the system until the
/// specified time.
///
/// # Example
///
/// ```
/// use rebop::define_system;
///
/// define_system! {
///     rtx rtl rdi rdm rdp;
///     Dimers { gene, mRNA, protein, dimer }
///     transcription   : gene      => gene + mRNA      @ rtx
///     translation     : mRNA      => mRNA + protein   @ rtl
///     dimerization    : 2 protein => dimer            @ rdi
///     decay_mRNA      : mRNA      =>                  @ rdm
///     decay_prot      : protein   =>                  @ rdp
/// }
/// let mut dimers = Dimers::with_parameters(25., 1000., 0.001, 0.1, 1.);
/// //                                       rtx,   rtl,   rdi, rdm,rdp
/// // but we could also have done
/// // let mut dimers = Dimers::new();
/// // dimers.rtx = 25.;
/// // etc.
/// dimers.gene = 1;
/// dimers.advance_until(1.);
/// println!("t = {}, dimer = {}", dimers.t, dimers.dimer);
/// ```
#[macro_export]
macro_rules! define_system {
    (
      $($param:ident)*;
      $name:ident { $($species:ident),* }
      $($rname:ident:
          $($($nr:literal)? $r:ident)? $(+ $($tnr:literal)? $tr:ident)* =>
          $($($np:literal)? $p:ident)? $(+ $($tnp:literal)? $tp:ident)*
          @ $rate:expr)*
      ) => {
        /// Structure representing the problem, with the species and the time.
        #[allow(non_snake_case)]
        #[derive(Clone, Debug)]
        struct $name {
            $($species:isize,)*
            $($param:f64,)*
            t: f64,
            rng: $crate::rand::rngs::SmallRng,
        }
        impl $name {
            /// Constructs an object representing the problem.
            fn new() -> Self {
                use $crate::rand::SeedableRng;
                $name {
                    $($species: 0,)*
                    $($param: f64::NAN,)*
                    t: 0.,
                    rng: $crate::rand::rngs::SmallRng::from_os_rng()
                }
            }
            /// Seeds the random number generator.
            fn seed(&mut self, seed: u64) {
                use $crate::rand::SeedableRng;
                self.rng = $crate::rand::rngs::SmallRng::seed_from_u64(seed);
            }
            /// Constructs an object representing the problem,
            /// specifying parameter values.
            #[allow(non_snake_case)]
            fn with_parameters($($param: f64),*) -> Self {
                use $crate::rand::SeedableRng;
                $name {
                    $($species: 0,)*
                    $($param,)*
                    t: 0.,
                    rng: $crate::rand::rngs::SmallRng::from_os_rng()
                }
            }
            /// Simulates the problem until `t = tmax`.
            #[allow(non_snake_case)]
            fn advance_until(&mut self, tmax: f64) {
                use $crate::rand::Rng;
                $(let $param = self.$param;)*
                $(let $species = self.$species as f64;)*
                loop {
                    $(let $rname = $rate $(* $crate::_rate_lma!($($nr)? * self.$r))? $(* $crate::_rate_lma!($($tnr)? * self.$tr) )*;)*
                    let total_rate = 0. $(+ $rname)*;
                    // we don't want to use partial_cmp, for performance
                    #[allow(clippy::neg_cmp_op_on_partial_ord)]
                    if !(total_rate > 0.) {
                        self.t = tmax;
                        return
                    }
                    self.t += self.rng.sample::<f64, _>($crate::rand_distr::Exp1) / total_rate;
                    if self.t > tmax {
                        self.t = tmax;
                        return
                    }
                    let reaction_choice = total_rate * self.rng.random::<f64>();
                    $crate::_choice!(self reaction_choice 0.;
                        $($rname:
                            $($($nr)? $r)? $(+ $($tnr)? $tr)* =>
                            $($($np)? $p)? $(+ $($tnp)? $tp)*;)*);
                }
            }
        }
    };
}

/// Auxiliary macro used in `define_system`.
#[macro_export]
macro_rules! _rate_lma {
    ($($n:literal)? * $species:expr) => {
        {
            let mut rate = $species;
            $(
                for i in 1..$n {
                    rate *= $species - i;
                }
            )?
            rate as f64
        }
    }
}

/// Auxiliary macro used in `define_system`.
#[macro_export]
macro_rules! _choice {
    ($self:ident $rc:ident $carry:expr; ) => {};
    ($self:ident $rc:ident $carry:expr;
     $rname:ident:
     $($($nr:literal)? $r:ident)? $(+ $($tnr:literal)? $tr:ident)* =>
     $($($np:literal)? $p:ident)? $(+ $($tnp:literal)? $tp:ident)*;
     $($tail:ident:
         $($($tailnr:literal)? $tailr:ident)? $(+ $($tailtnr:literal)? $tailtr:ident)* =>
         $($($tailnp:literal)? $tailp:ident)? $(+ $($tailtnp:literal)? $tailtp:ident)*;)*) => {
        if $rc < $carry + $rname {
            $($self.$r -= 1 $(+ $nr - 1)?;)?
            $($self.$tr -= 1 $(+ $tnr - 1)?;)*
            $($self.$p += 1 $(+ $np - 1)?;)?
            $($self.$tp += 1 $(+ $tnp - 1)?;)*
        } else {
            $crate::_choice!($self $rc $carry + $rname;
                $($tail:
                    $($($tailnr)? $tailr)? $(+ $($tailtnr)? $tailtr)* =>
                    $($($tailnp)? $tailp)? $(+ $($tailtnp)? $tailtp)*;)*);
        }
    };
}

#[cfg(test)]
mod tests {
    #[test]
    fn sir() {
        define_system! {
            r1 r2;
            SIR { S, I, R }
            r_infection: S + I  => I + I    @ r1
            r_remission: I      => R        @ r2
        }
        let mut sir = SIR::new();
        sir.r1 = 0.1 / 10000.;
        sir.r2 = 0.01;
        sir.S = 9999;
        sir.I = 1;
        sir.advance_until(1000.);
        assert_eq!(sir.S + sir.I + sir.R, 10000);
    }
    #[test]
    fn dimers() {
        define_system! {
            rtx rtl rdi rdm rdp;
            Dimers { gene, mRNA, protein, dimer }
            r_tx : gene         => gene + mRNA      @ rtx
            r_tl : mRNA         => mRNA + protein   @ rtl
            r_di : 2 protein    => dimer            @ rdi
            r_dm : mRNA         =>                  @ rdm
            r_dp : protein      =>                  @ rdp
        }
        let mut dimers = Dimers::with_parameters(25., 1000., 0.001, 0.1, 1.);
        dimers.gene = 1;
        dimers.advance_until(1.);
        assert_eq!(dimers.gene, 1);
        assert!(1000 < dimers.dimer);
        assert!(dimers.dimer < 10000);
    }
    #[test]
    fn birth_death() {
        define_system! {
            r_birth r_death;
            BirthDeath { A }
            birth:      => A    @ r_birth
            death:  A   =>      @ r_death
        }
        let mut birth_death = BirthDeath::new();
        birth_death.r_birth = 10.;
        birth_death.r_death = 0.1;
        birth_death.advance_until(100.);
        assert!(50 < birth_death.A);
        assert!(birth_death.A < 200);
    }
    #[test]
    fn birth_death_forgot_a_parameter() {
        define_system! {
            r_birth r_death;
            BirthDeath { A }
            birth:      => A    @ r_birth
            death:  A   =>      @ r_death
        }
        let mut birth_death = BirthDeath::new();
        birth_death.r_birth = 10.;
        // forgot to define r_death which is then NaN, so no reaction happens
        birth_death.advance_until(100.);
        assert!((birth_death.t - 100.).abs() < f64::EPSILON);
        assert_eq!(birth_death.A, 0);
    }
    #[test]
    fn no_reactions() {
        define_system! {
            ;
            FooBarBuz { Foo, Bar, Buz }
        }
        let mut foobarbuz = FooBarBuz::new();
        foobarbuz.Foo = 42;
        foobarbuz.Bar = 1337;
        foobarbuz.advance_until(1e20);
        assert!((foobarbuz.t - 1e20).abs() < f64::EPSILON);
        assert_eq!(foobarbuz.Foo, 42);
        assert_eq!(foobarbuz.Bar, 1337);
        assert_eq!(foobarbuz.Buz, 0);
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
