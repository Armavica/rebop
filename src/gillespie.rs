//! Function-based API to describe chemical reaction networks and
//! simulate them.

use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::Exp1;
use std::cmp::Ordering;

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
pub struct Gillespie<T: AsIndex, const N: usize> {
    species: [isize; N],
    t: f64,
    rates: Vec<Rate<T>>,
    reactions: Vec<[isize; N]>,
    rng: SmallRng,
}

impl<T: AsIndex + Clone, const N: usize> Gillespie<T, N> {
    /// Creates a new problem instance, with `N` different species of
    /// specified initial conditions.
    pub fn new(species: [isize; N]) -> Self {
        Gillespie {
            species,
            t: 0.,
            rates: Vec::new(),
            reactions: Vec::new(),
            rng: SmallRng::from_entropy(),
        }
    }
    /// Returns the number of species in the problem.
    ///
    /// ```
    /// use rebop::gillespie::Gillespie;
    /// let mut p: Gillespie<usize, 4> = Gillespie::new([0, 1, 10, 100]);
    /// assert_eq!(p.nb_species(), 4);
    /// ```
    pub fn nb_species(&self) -> usize {
        N
    }
    /// Returns the number of reactions in the problem.
    ///
    /// ```
    /// use rebop::gillespie::Gillespie;
    /// let mut p: Gillespie<usize, 4> = Gillespie::new([0, 1, 10, 100]);
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
    pub fn add_reaction(&mut self, rate: Rate<T>, reaction: [isize; N]) {
        self.rates.push(rate);
        self.reactions.push(reaction);
    }
    /// Returns the current time in the model.
    pub fn get_time(&self) -> f64 {
        self.t
    }
    /// Returns the current amount of a species.
    ///
    /// ```
    /// use rebop::gillespie::Gillespie;
    /// let p: Gillespie<usize, 4> = Gillespie::new([0, 1, 10, 100]);
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
    /// dimers.add_reaction(Rate::new(0.001, &[SRate::LMA(Dimers::P), SRate::LMA(Dimers::P)]), [0, 0, -2, 1]);
    /// dimers.add_reaction(Rate::new(0.1, &[SRate::LMA(Dimers::M)]), [0, -1, 0, 0]);
    /// dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::P)]), [0, 0, -1, 0]);
    /// assert_eq!(dimers.get_time(), 0.);
    /// assert_eq!(dimers.get_species(&Dimers::D), 0);
    /// dimers.advance_until(1.);
    /// assert_eq!(dimers.get_time(), 1.);
    /// assert!(dimers.get_species(&Dimers::D) > 0);
    /// ```
    pub fn advance_until(&mut self, tmax: f64) {
        let mut cumulative_rates = Vec::with_capacity(self.rates.len());
        while self.t < tmax {
            let mut total_rate = 0.;
            for r in &self.rates {
                total_rate += r.rate(&self.species);
                cumulative_rates.push(total_rate);
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
            let chosen_rate = total_rate * self.rng.gen::<f64>();
            let ireaction = match cumulative_rates.binary_search_by(|w| {
                if *w <= chosen_rate {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            }) {
                Ok(ireaction) | Err(ireaction) => ireaction,
            };
            // here we have ireaction < self.reactions.len() because chosen_rate < total_rate
            // FIXME: remove the bound check
            for (i, &r) in self.reactions[ireaction].iter().enumerate() {
                self.species[i] += r;
            }
            cumulative_rates.clear();
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
#[derive(Clone, Debug)]
pub enum SRate<T: AsIndex> {
    /// Law of mass action: `LMA(T) = [T]`
    LMA(T),
    /// Michaelis--Menten: `MM(T, KC) = [T] / (KC + [T])`
    MM(T, f64),
    /// Positive Hill function: `PosHill(T, KC, n) = [T]^n / (KC^n + [T]^n)`
    PosHill(T, f64, f64),
    /// Negative Hill function: `NegHill(T, KC, n) = 1 - [T]^n / (KC^n + [T]^n)`
    NegHill(T, f64, f64),
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
        dimers.add_reaction(
            Rate::new(0.001, &[SRate::LMA(Dimers::P), SRate::LMA(Dimers::P)]),
            [0, 0, -2, 1],
        );
        dimers.add_reaction(Rate::new(0.1, &[SRate::LMA(Dimers::M)]), [0, -1, 0, 0]);
        dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::P)]), [0, 0, -1, 0]);
        dimers.advance_until(1.);
        println!("{:?}", dimers.species);
        assert_eq!(dimers.get_species(&Dimers::G), 1);
        assert!(1000 < dimers.get_species(&Dimers::D));
        assert!(dimers.get_species(&Dimers::D) < 10000);
    }
}
