//! Function-based API to describe chemical reaction networks and
//! simulate them.

use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::Exp1;

#[derive(Clone, Debug)]
pub enum Expr {
    Constant(f64),
    Concentration(usize),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
    Pow(Box<Expr>, Box<Expr>),
    Exp(Box<Expr>),
}

impl Expr {
    fn eval(&self, species: &[isize]) -> f64 {
        match self {
            Expr::Constant(c) => *c,
            Expr::Concentration(i) => *unsafe { species.get_unchecked(*i) } as f64,
            Expr::Add(a, b) => a.eval(species) + b.eval(species),
            Expr::Sub(a, b) => a.eval(species) - b.eval(species),
            Expr::Mul(a, b) => a.eval(species) * b.eval(species),
            Expr::Div(a, b) => a.eval(species) / b.eval(species),
            Expr::Pow(a, b) => a.eval(species).powf(b.eval(species)),
            Expr::Exp(a) => a.eval(species).exp(),
        }
    }
}

#[derive(Clone, Debug)]
pub enum Rate {
    LMA(f64, Vec<u32>),
    LMASparse(f64, Vec<(u32, u32)>),
    Expr(Expr),
}

impl Rate {
    pub fn lma<V: AsRef<[u32]>>(rate: f64, reactants: V) -> Self {
        Rate::LMA(rate, reactants.as_ref().to_vec())
    }
    pub fn sparse(self) -> Self {
        match self {
            Rate::LMA(rate, reactants) => {
                let sparse = reactants
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &exponent)| {
                        Some((index as u32, exponent)).filter(|&(_, exponent)| exponent > 0)
                    })
                    .collect();
                Rate::LMASparse(rate, sparse)
            }
            Rate::LMASparse(_, _) => self,
            Rate::Expr(_) => unimplemented!(),
        }
    }
    fn rate(&self, species: &[isize]) -> f64 {
        match self {
            Rate::LMA(rate, ref reactants) => species
                .iter()
                .zip(reactants.iter())
                .fold(*rate, |acc, (&n, &e)| {
                    (n + 1 - e as isize..=n).fold(acc, |acc, x| acc * x as f64)
                }),
            Rate::LMASparse(mut rate, sparse) => {
                for &(index, exponent) in sparse.iter() {
                    let n = *unsafe { species.get_unchecked(index as usize) };
                    for i in (n + 1 - exponent as isize)..=n {
                        rate *= i as f64;
                    }
                }
                rate
            }
            Rate::Expr(expr) => expr.eval(species),
        }
    }
}

#[derive(Clone, Debug)]
pub enum Jump {
    Flat(Vec<isize>),
    Sparse(Vec<(usize, isize)>),
}

impl Jump {
    pub fn new<V: AsRef<[isize]>>(differences: V) -> Self {
        Jump::Flat(differences.as_ref().to_vec())
    }
    pub fn new_sparse<V: AsRef<[(usize, isize)]>>(sparse: V) -> Self {
        Jump::Sparse(sparse.as_ref().to_vec())
    }
    pub fn sparse(self) -> Self {
        match self {
            Jump::Flat(differences) => {
                let sparse = differences
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &difference)| {
                        Some((index, difference)).filter(|&(_, difference)| difference != 0)
                    })
                    .collect();
                Jump::Sparse(sparse)
            }
            Jump::Sparse(_) => self,
        }
    }
    fn affect(&self, species: &mut [isize]) {
        match self {
            Jump::Flat(differences) => species
                .iter_mut()
                .zip(differences.iter())
                .for_each(|(s, d)| *s += d),
            Jump::Sparse(differences) => differences.iter().for_each(|&(index, difference)| {
                *unsafe { species.get_unchecked_mut(index) } += difference
            }),
        }
    }
}

/// Main structure, represents the problem and contains simulation methods.
#[derive(Clone, Debug)]
pub struct Gillespie {
    species: Vec<isize>,
    t: f64,
    reactions: Vec<(Rate, Jump)>,
    rng: SmallRng,
}

impl Gillespie {
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
    pub fn new_with_seed<V: AsRef<[isize]>>(species: V, seed: u64) -> Self {
        Gillespie {
            species: species.as_ref().to_vec(),
            t: 0.,
            reactions: Vec::new(),
            rng: SmallRng::seed_from_u64(seed),
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
    /// let mut p: Gillespie = Gillespie::new([0, 1, 10, 100]);
    /// assert_eq!(p.nb_species(), 4);
    /// ```
    pub fn nb_species(&self) -> usize {
        self.species.len()
    }
    /// Returns the number of reactions in the problem.
    ///
    /// ```
    /// use rebop::gillespie::Gillespie;
    /// let mut p: Gillespie = Gillespie::new([0, 1, 10, 100]);
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
    /// use rebop::gillespie::{Gillespie, Rate};
    /// let mut sir = Gillespie::new([9999, 1, 0]);
    /// //                           [   S, I, R]
    /// // S + I -> I + I with rate 1e-5
    /// sir.add_reaction(Rate::lma(1e-5, [1, 1, 0]), [-1, 1, 0]);
    /// // I -> R with rate 0.01
    /// sir.add_reaction(Rate::lma(0.01, [0, 1, 0]), [0, -1, 1]);
    /// ```
    pub fn add_reaction<V: AsRef<[isize]>>(&mut self, rate: Rate, differences: V) {
        // This assert ensures that the jump does not go out of bounds of the species
        assert_eq!(differences.as_ref().len(), self.species.len());
        let jump = Jump::new(differences);
        self.reactions.push((rate.sparse(), jump));
    }
    /// Returns the current time in the model.
    pub fn get_time(&self) -> f64 {
        self.t
    }
    /// Sets the current time in the model.
    pub fn set_time(&mut self, t: f64) {
        self.t = t;
    }
    /// Returns the current amount of a species.
    ///
    /// ```
    /// use rebop::gillespie::Gillespie;
    /// let p: Gillespie = Gillespie::new([0, 1, 10, 100]);
    /// assert_eq!(p.get_species(2), 10);
    /// ```
    pub fn get_species(&self, s: usize) -> isize {
        self.species[s]
    }
    /// Sets the amount of species in the model.
    pub fn set_species<V: AsRef<[isize]>>(&mut self, species: V) {
        assert_eq!(species.as_ref().len(), self.species.len());
        self.species = species.as_ref().to_vec();
    }
    /// Simulates the problem until `tmax`.
    ///
    /// ```
    /// use rebop::gillespie::{Gillespie, Rate};
    /// let mut dimers = Gillespie::new([1, 0, 0, 0]);
    /// //                              [G, M, P, D]
    /// dimers.add_reaction(Rate::lma(25., [1, 0, 0, 0]), [0, 1, 0, 0]);
    /// dimers.add_reaction(Rate::lma(1000., [0, 1, 0, 0]), [0, 0, 1, 0]);
    /// dimers.add_reaction(Rate::lma(0.001, [0, 0, 2, 0]), [0, 0, -2, 1]);
    /// dimers.add_reaction(Rate::lma(0.1, [0, 1, 0, 0]), [0, -1, 0, 0]);
    /// dimers.add_reaction(Rate::lma(1., [0, 0, 1, 0]), [0, 0, -1, 0]);
    /// assert_eq!(dimers.get_time(), 0.);
    /// assert_eq!(dimers.get_species(3), 0);
    /// dimers.advance_until(1.);
    /// assert_eq!(dimers.get_time(), 1.);
    /// assert!(dimers.get_species(3) > 0);
    /// ```
    pub fn advance_until(&mut self, tmax: f64) {
        let mut rates = vec![f64::NAN; self.reactions.len()];
        loop {
            //let total_rate = make_rates(&self.reactions, &self.species, &mut rates);
            let total_rate = make_cumrates(&self.reactions, &self.species, &mut rates);

            // we don't want to use partial_cmp, for performance
            #[allow(clippy::neg_cmp_op_on_partial_ord)]
            if !(0. < total_rate) {
                self.t = tmax;
                return;
            }
            self.t += self.rng.sample::<f64, _>(Exp1) / total_rate;
            if self.t > tmax {
                self.t = tmax;
                return;
            }
            let chosen_rate = total_rate * self.rng.gen::<f64>();

            //let ireaction = choose_rate_sum(chosen_rate, &rates);
            //let ireaction = choose_rate_for(chosen_rate, &rates);
            let ireaction = choose_cumrate_sum(chosen_rate, &rates);
            //let ireaction = choose_cumrate_for(chosen_rate, &rates);
            //let ireaction = choose_cumrate_takewhile(chosen_rate, &rates);
            // here we have ireaction < self.reactions.len() because chosen_rate < total_rate
            let reaction = unsafe { self.reactions.get_unchecked(ireaction) };

            reaction.1.affect(&mut self.species);
        }
    }
}

fn make_rates(reactions: &[(Rate, Jump)], species: &[isize], rates: &mut [f64]) -> f64 {
    let mut total_rate = 0.0;
    for ((rate, _), num_rate) in reactions.iter().zip(rates.iter_mut()) {
        *num_rate = rate.rate(species);
        total_rate += *num_rate;
    }
    total_rate
}

fn make_cumrates(reactions: &[(Rate, Jump)], species: &[isize], cum_rates: &mut [f64]) -> f64 {
    let mut total_rate = 0.0;
    for ((rate, _), cum_rate) in reactions.iter().zip(cum_rates.iter_mut()) {
        *cum_rate = total_rate + rate.rate(species);
        total_rate = *cum_rate;
    }
    total_rate
}

fn choose_rate_for(mut chosen_rate: f64, rates: &[f64]) -> usize {
    let mut ireaction = rates.len() - 1;
    for (ir, &rate) in rates.iter().enumerate() {
        chosen_rate -= rate;
        if chosen_rate < 0. {
            ireaction = ir;
            break;
        }
    }
    ireaction
}

fn choose_cumrate_for(chosen_rate: f64, cumrates: &[f64]) -> usize {
    let mut ireaction = cumrates.len() - 1;
    for (ir, &cumrate) in cumrates.iter().enumerate() {
        if chosen_rate < cumrate {
            ireaction = ir;
            break;
        }
    }
    ireaction
}

fn choose_rate_sum(chosen_rate: f64, rates: &[f64]) -> usize {
    rates
        .iter()
        .scan(0.0, |cum, &r| {
            *cum += r;
            Some(if *cum < chosen_rate { 1 } else { 0 })
        })
        .sum()
}

fn choose_cumrate_sum(chosen_rate: f64, cumrates: &[f64]) -> usize {
    cumrates
        .iter()
        .map(|&cum| if cum < chosen_rate { 1 } else { 0 })
        .sum()
}

fn choose_cumrate_takewhile(chosen_rate: f64, cumrates: &[f64]) -> usize {
    cumrates
        .iter()
        .take_while(|&&cum| cum < chosen_rate)
        .count()
}

#[cfg(test)]
mod tests {
    use crate::gillespie::{Gillespie, Rate};
    #[test]
    fn sir() {
        let mut sir = Gillespie::new([9999, 1, 0]);
        sir.add_reaction(Rate::lma(0.1 / 10000., [1, 1, 0]), [-1, 1, 0]);
        sir.add_reaction(Rate::lma(0.01, [0, 1, 0]), [0, -1, 1]);
        sir.advance_until(250.);
        assert_eq!(
            sir.get_species(0) + sir.get_species(1) + sir.get_species(2),
            10000
        );
    }
    #[test]
    fn dimers() {
        let mut dimers = Gillespie::new([1, 0, 0, 0]);
        dimers.add_reaction(Rate::lma(25., [1, 0, 0, 0]), [0, 1, 0, 0]);
        dimers.add_reaction(Rate::lma(1000., [0, 1, 0, 0]), [0, 0, 1, 0]);
        dimers.add_reaction(Rate::lma(0.001, [0, 0, 2, 0]), [0, 0, -2, 1]);
        dimers.add_reaction(Rate::lma(0.1, [0, 1, 0, 0]), [0, -1, 0, 0]);
        dimers.add_reaction(Rate::lma(1., [0, 0, 1, 0]), [0, 0, -1, 0]);
        dimers.advance_until(1.);
        assert_eq!(dimers.get_species(0), 1);
        assert!(1000 < dimers.get_species(2));
        assert!(dimers.get_species(3) < 10000);
    }
}
