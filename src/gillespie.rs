use rand::distributions::Distribution;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::Exp;
use std::cmp::Ordering;

pub trait AsIndex {
    fn as_index(&self) -> usize;
}

#[macro_export]
macro_rules! index_enum {
    (enum $name:ident $body:tt) => {
        #[allow(dead_code)]
        #[allow(non_camel_case_types)]
        #[derive(Copy, Clone)] enum $name $body
        impl AsIndex for $name { fn as_index(&self) -> usize { *self as usize } }
    }
}

pub struct Gillespie<T: AsIndex, const N: usize> {
    species: [isize; N],
    t: f64,
    rates: Vec<Rate<T>>,
    reactions: Vec<[isize; N]>,
    rng: SmallRng,
}

impl<T: AsIndex + Clone, const N: usize> Gillespie<T, N> {
    pub fn new(species: [isize; N]) -> Self {
        Gillespie {
            species,
            t: 0.,
            rates: Vec::new(),
            reactions: Vec::new(),
            rng: SmallRng::from_entropy(),
        }
    }
    pub fn nb_species(&self) -> usize {
        self.species.len()
    }
    pub fn nb_reactions(&self) -> usize {
        self.reactions.len()
    }
    pub fn add_reaction(&mut self, rate: Rate<T>, reaction: [isize; N]) {
        self.rates.push(rate);
        self.reactions.push(reaction);
    }
    pub fn get_species(&self, s: &T) -> isize {
        self.species[s.as_index()]
    }
    pub fn advance_until(&mut self, tmax: f64) {
        while self.t < tmax {
            let mut total_rate = 0.;
            let mut cumulative_rates = Vec::with_capacity(self.rates.len());
            for r in &self.rates {
                total_rate += r.rate(&self.species);
                cumulative_rates.push(total_rate);
            }
            // we don't want to use partial_cmp, for performance
            #[allow(clippy::neg_cmp_op_on_partial_ord)]
            if !(total_rate > 0.) {
                self.t = tmax;
                return
            }
            // unwrap is safe because we just checked that total_rate > 0
            self.t += Exp::new(total_rate).unwrap().sample(&mut self.rng);
            if self.t > tmax {
                self.t = tmax;
                return
            }
            let chosen_rate = total_rate * self.rng.gen::<f64>();
            let ireaction = match cumulative_rates.binary_search_by(
                |w| if *w <= chosen_rate { Ordering::Less } else { Ordering::Greater }) {
                Ok(ireaction) | Err(ireaction) => ireaction,
            };
            // here we have ireaction < self.reactions.len() because chosen_rate < total_rate
            // FIXME: remove the bound check
            for (i, &r) in self.reactions[ireaction].iter().enumerate() {
                self.species[i] += r;
            }
        }
    }
}

pub struct Rate<T: AsIndex> {
    rate: f64,
    species: Vec<SRate<T>>,
}

#[derive(Clone)]
pub enum SRate<T: AsIndex> {
    LMA(T),               // S
    LMA2(T),              // S
    MM(T, f64),           // S, KC
    PosHill(T, f64, f64), // S, KC, n
    NegHill(T, f64, f64), // S, KC, n
}

impl<T: AsIndex + Clone> Rate<T> {
    pub fn new(rate: f64, species: &[SRate<T>]) -> Self {
        Rate {
            rate,
            species: species.to_vec(),
        }
    }
    pub fn rate(&self, species: &[isize]) -> f64 {
        let mut r = self.rate;
        for factor in &self.species {
            match *factor {
                SRate::LMA(ref s) => r *= species[s.as_index()] as f64,
                SRate::LMA2(ref s) => r *= (species[s.as_index()] * species[s.as_index()]) as f64,
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
