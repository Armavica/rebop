//! Function-based API to describe chemical reaction networks and
//! simulate them.

use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::Exp1;

/// Trait that a species type must implement.
///
/// It makes sure that the species type can be used as an index.
pub trait AsIndex {
    fn as_index(&self) -> usize;
}

impl AsIndex for usize {
    fn as_index(&self) -> usize {
        *self
    }
}

/// Utility macro to automatically implement the [`AsIndex`] trait on an
/// enum.
#[macro_export]
macro_rules! index_enum {
    (enum $name:ident $body:tt) => {
        #[allow(dead_code)]
        #[allow(non_camel_case_types)]
        #[derive(Copy, Clone, Debug)] enum $name $body
        impl AsIndex for $name { fn as_index(&self) -> usize { *self as usize } }
    }
}

/// Main structure, represents the problem and contains simulation methods.
#[derive(Clone, Debug)]
pub struct Gillespie<T: AsIndex> {
    species: Vec<isize>,
    t: f64,
    reactions: Vec<(Rate<T>, Vec<isize>)>,
    rng: SmallRng,
}

impl<T: AsIndex + Clone> Gillespie<T> {
    /// Creates a new problem instance, with `N` different species of
    /// specified initial conditions.
    pub fn new<V: AsRef<[isize]>>(species: V) -> Self {
        Gillespie {
            species: species.as_ref().to_vec(),
            t: 0.,
            reactions: Vec::new(),
            rng: SmallRng::from_entropy(),
        }
    }
    /// Seeds the random number generator.
    pub fn seed(&mut self, seed: u64) {
        self.rng = SmallRng::seed_from_u64(seed);
    }
    /// Returns the number of species in the problem.
    ///
    /// ```
    /// use rebop::gillespie::Gillespie;
    /// let mut p: Gillespie<usize> = Gillespie::new([0, 1, 10, 100]);
    /// assert_eq!(p.nb_species(), 4);
    /// ```
    pub fn nb_species(&self) -> usize {
        self.species.len()
    }
    /// Returns the number of reactions in the problem.
    ///
    /// ```
    /// use rebop::gillespie::Gillespie;
    /// let mut p: Gillespie<usize> = Gillespie::new([0, 1, 10, 100]);
    /// assert_eq!(p.nb_reactions(), 0);
    /// ```
    pub fn nb_reactions(&self) -> usize {
        self.reactions.len()
    }
    /// Adds a reaction to the problem.
    ///
    /// `rate` is the reaction rate and `reaction` is an array
    /// describing the state change as a result of the reaction.
    /// ```
    /// use rebop::index_enum;
    /// use rebop::gillespie::{AsIndex, Gillespie, Rate, SRate};
    /// index_enum! { enum SIR { S, I, R } }
    /// let mut sir = Gillespie::new([9999, 1, 0]);
    /// // S + I -> I + I with rate 1e-5
    /// sir.add_reaction(
    ///     Rate::new(1e-5, &[SRate::LMA(SIR::S), SRate::LMA(SIR::I)]),
    ///     [-1, 1, 0],
    /// );
    /// // I -> R with rate 0.01
    /// sir.add_reaction(Rate::new(0.01, &[SRate::LMA(SIR::I)]), [0, -1, 1]);
    /// ```
    pub fn add_reaction<V: AsRef<[isize]>>(&mut self, rate: Rate<T>, reaction: V) {
        assert_eq!(reaction.as_ref().len(), self.species.len());
        self.reactions.push((rate, reaction.as_ref().to_vec()));
    }
    /// Returns the current time in the model.
    pub fn get_time(&self) -> f64 {
        self.t
    }
    /// Returns the current amount of a species.
    ///
    /// ```
    /// use rebop::gillespie::Gillespie;
    /// let p: Gillespie<usize> = Gillespie::new([0, 1, 10, 100]);
    /// assert_eq!(p.get_species(&2), 10);
    /// ```
    pub fn get_species(&self, s: &T) -> isize {
        self.species[s.as_index()]
    }
    /// Simulates the problem until `tmax`.
    ///
    /// ```
    /// use rebop::index_enum;
    /// use rebop::gillespie::{AsIndex, Gillespie, Rate, SRate};
    /// index_enum! { enum Dimers { G, M, P, D } }
    /// let mut dimers = Gillespie::new([1, 0, 0, 0]);
    /// dimers.add_reaction(Rate::new(25., &[SRate::LMA(Dimers::G)]), [0, 1, 0, 0]);
    /// dimers.add_reaction(Rate::new(1000., &[SRate::LMA(Dimers::M)]), [0, 0, 1, 0]);
    /// dimers.add_reaction(Rate::new(0.001, &[SRate::LMA2(Dimers::P)]), [0, 0, -2, 1]);
    /// dimers.add_reaction(Rate::new(0.1, &[SRate::LMA(Dimers::M)]), [0, -1, 0, 0]);
    /// dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::P)]), [0, 0, -1, 0]);
    /// assert_eq!(dimers.get_time(), 0.);
    /// assert_eq!(dimers.get_species(&Dimers::D), 0);
    /// dimers.advance_until(1.);
    /// assert_eq!(dimers.get_time(), 1.);
    /// assert!(dimers.get_species(&Dimers::D) > 0);
    /// ```
    pub fn advance_until(&mut self, tmax: f64) {
        let mut rates = vec![f64::NAN; self.reactions.len()];
        loop {
            let mut total_rate = 0.;
            for ((rate, _), num_rate) in self.reactions.iter().zip(rates.iter_mut()) {
                *num_rate = rate.rate(&self.species);
                total_rate += *num_rate;
            }
            // we don't want to use partial_cmp, for performance
            #[allow(clippy::neg_cmp_op_on_partial_ord)]
            if !(total_rate > 0.) {
                self.t = tmax;
                return;
            }
            self.t += self.rng.sample::<f64, _>(Exp1) / total_rate;
            if self.t > tmax {
                self.t = tmax;
                return;
            }
            let mut chosen_rate = total_rate * self.rng.gen::<f64>();
            let mut ireaction = self.reactions.len() - 1;
            for (ir, rate) in rates.iter().enumerate() {
                chosen_rate -= rate;
                if chosen_rate < 0. {
                    ireaction = ir;
                    break
                }
            };
            // here we have ireaction < self.reactions.len() because chosen_rate < total_rate
            // FIXME: remove the bound check
            for (i, &r) in self.reactions[ireaction].1.iter().enumerate() {
                self.species[i] += r;
            }
        }
    }
}

/// This structure represents a reaction rate, as the product of
/// a numerical constant and several factors, each represented by
/// a [`SRate`].
#[derive(Clone, Debug)]
pub struct Rate<T: AsIndex> {
    rate: f64,
    species: Vec<SRate<T>>,
}

/// This enum represents a factor of a reaction rate, used in the
/// construction of a [`Rate`].
#[derive(Clone, Debug, PartialEq)]
pub enum SRate<T: AsIndex> {
    /// Law of mass action: `LMA(T) = [T]`
    LMA(T),
    /// Law of mass action of order 2: `LMA2(T) = [T] * ([T] - 1)`
    LMA2(T),
    /// Law of mass action of order n: `LMAn(T, n) = [T] * ([T] - 1) * ... * ([T] - n + 1)`
    LMAn(T, usize),
    /// Michaelis--Menten: `MM(T, KC) = [T] / (KC + [T])`
    MM(T, f64),
    /// Positive Hill function: `PosHill(T, KC, n) = [T]^n / (KC^n + [T]^n)`
    PosHill(T, f64, f64),
    /// Negative Hill function: `NegHill(T, KC, n) = 1 - [T]^n / (KC^n + [T]^n)`
    NegHill(T, f64, f64),
}

impl<T: AsIndex> SRate<T> {
    fn as_index(&self) -> usize {
        match self {
            SRate::LMA(x) | SRate::LMA2(x) | SRate::LMAn(x, _) | SRate::MM(x, _) | SRate::PosHill(x, _, _) | SRate::NegHill(x, _, _) => x.as_index(),
        }
    }
}

impl<T: AsIndex + Clone> Rate<T> {
    /// Creates a reaction rate with a numerical constant and several factors.
    ///
    /// ```
    /// use rebop::gillespie::{Rate, SRate};
    /// // r1 = 25 C₀ C₁
    /// let r1 = Rate::new(25., &[SRate::LMA(0), SRate::LMA(1)]);
    /// // r2 = 4 C₀ C₁ / (10 + C₀)
    /// let r2 = Rate::new(4., &[SRate::MM(0, 10.), SRate::LMA(1)]);
    /// ```
    pub fn new(rate: f64, species: &[SRate<T>]) -> Self {
        Rate {
            rate,
            species: species.to_vec(),
        }
    }
    fn rate(&self, species: &[isize]) -> f64 {
        let mut r = self.rate;
        for factor in &self.species {
            match *factor {
                SRate::LMA(ref s) => r *= species[s.as_index()] as f64,
                SRate::LMA2(ref s) => r *= species[s.as_index()] as f64 * (species[s.as_index()] - 1) as f64,
                SRate::LMAn(ref s, n) => for i in 0..n {
                    r *= (species[s.as_index()] - i as isize) as f64;
                },
                SRate::MM(ref s, k) => {
                    r *= species[s.as_index()] as f64 / (k + species[s.as_index()] as f64)
                }
                SRate::PosHill(ref s, k, n) => {
                    r *= (species[s.as_index()] as f64).powf(n)
                        / (k.powf(n) + (species[s.as_index()] as f64).powf(n))
                }
                SRate::NegHill(ref s, k, n) => {
                    r *= (1. + (species[s.as_index()] as f64 / k).powf(n)).recip()
                }
            }
        }
        r
    }
}

#[cfg(test)]
mod tests {
    use crate::gillespie::{AsIndex, Gillespie, Rate, SRate};
    #[test]
    fn sir() {
        index_enum! { enum SIR { S, I, R } }
        let mut sir = Gillespie::new([9999, 1, 0]);
        sir.add_reaction(
            Rate::new(0.1 / 10000., &[SRate::LMA(SIR::S), SRate::LMA(SIR::I)]),
            [-1, 1, 0],
        );
        sir.add_reaction(Rate::new(0.01, &[SRate::LMA(SIR::I)]), [0, -1, 1]);
        sir.advance_until(250.);
        assert_eq!(
            sir.get_species(&SIR::S) + sir.get_species(&SIR::I) + sir.get_species(&SIR::R),
            10000
        );
    }
    #[test]
    fn dimers() {
        index_enum! { enum Dimers { G, M, P, D } }
        let mut dimers = Gillespie::new([1, 0, 0, 0]);
        dimers.add_reaction(Rate::new(25., &[SRate::LMA(Dimers::G)]), [0, 1, 0, 0]);
        dimers.add_reaction(Rate::new(1000., &[SRate::LMA(Dimers::M)]), [0, 0, 1, 0]);
        dimers.add_reaction(Rate::new(0.001, &[SRate::LMA2(Dimers::P)]), [0, 0, -2, 1]);
        dimers.add_reaction(Rate::new(0.1, &[SRate::LMA(Dimers::M)]), [0, -1, 0, 0]);
        dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::P)]), [0, 0, -1, 0]);
        dimers.advance_until(1.);
        println!("{:?}", dimers.species);
        assert_eq!(dimers.get_species(&Dimers::G), 1);
        assert!(1000 < dimers.get_species(&Dimers::D));
        assert!(dimers.get_species(&Dimers::D) < 10000);
    }
}
