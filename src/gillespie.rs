//! Function-based API to describe chemical reaction networks and
//! simulate them.

use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::Exp1;

use crate::expr::Expr;

#[derive(Clone, Debug, PartialEq)]
pub enum Rate {
    LMA(f64, Vec<u32>),
    LMASparse(f64, Vec<(u32, u32)>),
    Expr(Expr),
}

impl Rate {
    pub fn lma<V: AsRef<[u32]>>(rate: f64, stoechiometries: V) -> Self {
        Rate::LMA(rate, stoechiometries.as_ref().to_vec())
    }
    pub fn lma_sparse<V: AsRef<[(u32, u32)]>>(rate: f64, stoechiometries: V) -> Self {
        Rate::LMASparse(rate, stoechiometries.as_ref().to_vec())
    }
    pub fn expr(expr: Expr) -> Self {
        Rate::Expr(expr)
    }
    pub fn dense(self) -> Self {
        match self {
            Rate::LMASparse(rate, stoechiometries) => {
                let n = stoechiometries
                    .iter()
                    .map(|(i, _)| i + 1)
                    .max()
                    .unwrap_or(0);
                let mut dense = vec![0; n as usize];
                for (index, exponent) in stoechiometries {
                    dense[index as usize] += exponent;
                }
                Rate::LMA(rate, dense)
            }
            Rate::LMA(_, _) | Rate::Expr(_) => self,
        }
    }
    pub fn sparse(self) -> Self {
        match self {
            Rate::LMA(rate, stoechiometries) => {
                let sparse = stoechiometries
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &exponent)| {
                        Some((index as u32, exponent)).filter(|&(_, exponent)| exponent > 0)
                    })
                    .collect();
                Rate::LMASparse(rate, sparse)
            }
            Rate::LMASparse(_, _) | Rate::Expr(_) => self,
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

#[derive(Clone, Debug, PartialEq)]
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
    pub fn dense(self) -> Self {
        match self {
            Jump::Sparse(sparse) => {
                let n = sparse.iter().map(|(i, _)| i + 1).max().unwrap_or(0);
                let mut dense = vec![0; n];
                for (index, diff) in sparse {
                    dense[index] += diff;
                }
                Jump::Flat(dense)
            }
            Jump::Flat(_) => self,
        }
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
    sparse: bool,
}

impl Gillespie {
    /// Creates a new problem instance, with `N` different species of
    /// specified initial conditions.
    pub fn new<V: AsRef<[isize]>>(species: V, sparse: bool) -> Self {
        Gillespie {
            species: species.as_ref().to_vec(),
            t: 0.,
            reactions: Vec::new(),
            rng: SmallRng::from_os_rng(),
            sparse,
        }
    }
    pub fn new_with_seed<V: AsRef<[isize]>>(species: V, sparse: bool, seed: u64) -> Self {
        Gillespie {
            species: species.as_ref().to_vec(),
            t: 0.,
            reactions: Vec::new(),
            rng: SmallRng::seed_from_u64(seed),
            sparse,
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
    /// let mut p: Gillespie = Gillespie::new([0, 1, 10, 100], false);
    /// assert_eq!(p.nb_species(), 4);
    /// ```
    pub fn nb_species(&self) -> usize {
        self.species.len()
    }
    /// Returns the number of reactions in the problem.
    ///
    /// ```
    /// use rebop::gillespie::Gillespie;
    /// let mut p: Gillespie = Gillespie::new([0, 1, 10, 100], false);
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
    /// let mut sir = Gillespie::new([9999, 1, 0], false);
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
        if self.sparse {
            self.reactions.push((rate.sparse(), jump.sparse()));
        } else {
            self.reactions.push((rate, jump));
        }
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
    /// let p: Gillespie = Gillespie::new([0, 1, 10, 100], false);
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
    /// Simulates the problem until the next discrete reaction.
    pub fn advance_one_reaction(&mut self) {
        let mut rates = vec![f64::NAN; self.nb_reactions()];
        self._advance_one_reaction(&mut rates);
    }

    #[inline]
    pub fn _advance_one_reaction(&mut self, rates: &mut [f64]) {
        // let total_rate = make_rates(&self.reactions, &self.species, rates);
        let total_rate = make_cumrates(&self.reactions, &self.species, rates);

        // we don't want to use partial_cmp, for performance
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(0. < total_rate) {
            self.t = f64::INFINITY;
            return;
        }
        self.t += self.rng.sample::<f64, _>(Exp1) / total_rate;
        let chosen_rate = total_rate * self.rng.random::<f64>();

        // let ireaction = choose_rate_sum(chosen_rate, &rates);
        // let ireaction = choose_rate_for(chosen_rate, &rates);
        let ireaction = choose_cumrate_sum(chosen_rate, &rates);
        // let ireaction = choose_cumrate_for(chosen_rate, &rates);
        // let ireaction = choose_cumrate_takewhile(chosen_rate, &rates);
        // here we have ireaction < self.reactions.len() because chosen_rate < total_rate
        let reaction = unsafe { self.reactions.get_unchecked(ireaction) };

        reaction.1.affect(&mut self.species);
    }
    /// Simulates the problem until `tmax`.
    ///
    /// ```
    /// use rebop::gillespie::{Gillespie, Rate};
    /// let mut dimers = Gillespie::new([1, 0, 0, 0], false);
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
            let chosen_rate = total_rate * self.rng.random::<f64>();

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
    use crate::gillespie::{Expr, Gillespie, Jump, Rate};

    #[test]
    fn rate_conversion() {
        let rate = Rate::lma(2.0, [3]);
        assert_eq!(rate.clone().dense(), rate);
        assert_eq!(rate.clone().sparse(), Rate::lma_sparse(2.0, [(0, 3)]));
        assert_eq!(rate.clone().sparse().dense(), rate);
        let rate = Rate::lma_sparse(2.0, [(0, 1), (3, 2)]);
        assert_eq!(rate.clone().sparse(), rate);
        assert_eq!(rate.clone().dense(), Rate::lma(2.0, [1, 0, 0, 2]));
        assert_eq!(rate.clone().dense().sparse(), rate);
        let rate = Rate::expr(Expr::Constant(4.1));
        assert_eq!(rate.clone().dense(), rate);
        assert_eq!(rate.clone().sparse(), rate);
    }

    #[test]
    fn rate_sparse() {
        let rate = Rate::lma(2.0, [3]);
        assert_eq!(rate.sparse(), Rate::LMASparse(2.0, vec![(0, 3)]));
        let rate = Rate::lma(2.0, [0, 3]);
        assert_eq!(rate.sparse(), Rate::LMASparse(2.0, vec![(1, 3)]));
        let rate = Rate::lma(2.0, [0, 3, 0]);
        assert_eq!(rate.sparse(), Rate::LMASparse(2.0, vec![(1, 3)]));
        let rate = Rate::lma(2.0, [0, 3, 0, 1, 0, 0]);
        assert_eq!(rate.sparse(), Rate::LMASparse(2.0, vec![(1, 3), (3, 1)]));
    }

    #[test]
    fn rate_lma() {
        fn check_rate(rate: &Rate, species: &[isize], expected: f64) {
            assert_eq!(rate.clone().dense().rate(species), expected);
            assert_eq!(rate.clone().sparse().rate(species), expected);
        }

        let species = [5, 3];
        check_rate(&Rate::lma(2.0, [1, 0]), &species, 10.0);
        check_rate(&Rate::lma(2.0, [2, 0]), &species, 40.0);
        check_rate(&Rate::lma(2.0, [3, 0]), &species, 120.0);
        check_rate(&Rate::lma(2.0, [4, 0]), &species, 240.0);
        check_rate(&Rate::lma(2.0, [5, 0]), &species, 240.0);
        check_rate(&Rate::lma(2.0, [6, 0]), &species, 0.0);
        check_rate(&Rate::lma(2.0, [7, 0]), &species, 0.0);
        check_rate(&Rate::lma(2.0, [0, 1]), &species, 6.0);
        check_rate(&Rate::lma(2.0, [0, 2]), &species, 12.0);
        check_rate(&Rate::lma(2.0, [0, 3]), &species, 12.0);
        check_rate(&Rate::lma(2.0, [0, 4]), &species, 0.0);
        check_rate(&Rate::lma(2.0, [0, 5]), &species, 0.0);
        check_rate(&Rate::lma(2.0, [1, 1]), &species, 30.0);
        check_rate(&Rate::lma(2.0, [1, 2]), &species, 60.0);
        check_rate(&Rate::lma(2.0, [2, 1]), &species, 120.0);
        check_rate(&Rate::lma(2.0, [1, 20]), &species, 0.0);
    }

    #[test]
    fn jump_conversion() {
        let djump = Jump::new([-1, 0, 3, -2]);
        let sjump = Jump::new_sparse([(0, -1), (2, 3), (3, -2)]);
        assert_eq!(djump.clone().dense(), djump);
        assert_eq!(djump.clone().sparse(), sjump);
        assert_eq!(sjump.clone().sparse(), sjump);
        assert_eq!(sjump.clone().dense(), djump);
    }

    #[test]
    fn jump_affect() {
        let djump = Jump::new([-1, 0, 3, -2]);
        let mut species = vec![10; 4];
        djump.affect(&mut species);
        assert_eq!(species, vec![9, 10, 13, 8]);
        let sjump = Jump::new_sparse([(0, -1), (2, 3), (3, -2)]);
        let mut species = vec![10; 4];
        sjump.affect(&mut species);
        assert_eq!(species, vec![9, 10, 13, 8]);
    }

    #[test]
    fn gillespie_sir() {
        let mut sir = Gillespie::new([9999, 1, 0], false);
        sir.add_reaction(Rate::lma(0.1 / 10000., [1, 1, 0]), [-1, 1, 0]);
        sir.add_reaction(Rate::lma(0.01, [0, 1, 0]), [0, -1, 1]);
        for i in 1..=250 {
            sir.advance_until(i as f64);
            assert_eq!(
                sir.get_species(0) + sir.get_species(1) + sir.get_species(2),
                10000
            );
        }
    }

    #[test]
    fn gillespie_dimers() {
        let mut dimers = Gillespie::new([1, 0, 0, 0], false);
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
