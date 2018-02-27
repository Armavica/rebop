use rand::distributions::{Exp, IndependentSample, Range};
use rand::thread_rng;

pub trait AsIndex {
    fn as_index(&self) -> usize;
}

#[macro_export]
macro_rules! index_enum {
    (enum $name:ident $body:tt) => {
        #[derive(Copy, Clone)] enum $name $body
        impl AsIndex for $name { fn as_index(&self) -> usize { *self as usize } }
    }
}

impl<T: AsIndex+Clone> Gillespie<T> {
    pub fn new(species: &[isize]) -> Self {
        Gillespie {
            species: species.to_vec(),
            t: 0.,
            rates: Vec::new(),
            reactions: Vec::new()
        }
    }
    pub fn nb_species(&self) -> usize {
        self.species.len()
    }
    pub fn nb_reactions(&self) -> usize {
        self.reactions.len()
    }
    pub fn add_reaction(&mut self, rate: Rate<T>, reaction: &[isize]) {
        self.rates.push(rate);
        self.reactions.push(reaction.to_vec());
    }
    pub fn get_species(&self, s: &T) -> isize {
        self.species[s.as_index()]
    }
    fn choose_reaction(&self) -> Option<(usize, f64)> {
        let rates: Vec<f64> = self.rates.iter().map(|r| r.rate(&self.species)).collect();
        let total_rate: f64 = rates.iter().sum();
        if total_rate <= 0. {
            return None
        }
        let mut rng = thread_rng();
        let exp = Exp::new(total_rate);
        let reaction_time = exp.ind_sample(&mut rng);
        let mut reaction_choice = Range::new(0., total_rate).ind_sample(&mut rng);
        for (i, &rate) in rates.iter().enumerate() {
            if rate >= reaction_choice {
                return Some((i, reaction_time))
            }
            reaction_choice -= rate;
        }
        unreachable!()
    }
    pub fn step(&mut self) -> bool {
        match self.choose_reaction() {
            Some((reaction, reaction_time)) => {
                for (i, &r) in self.reactions[reaction].iter().enumerate() {
                    self.species[i] += r;
                }
                self.t += reaction_time;
                true
            }
            None => false
        }
    }
    pub fn advance_until(&mut self, tmax: f64) -> Vec<(f64, Vec<isize>)> {
        let mut r = Vec::new();
        r.push((self.t, self.species.clone()));
        while self.t < tmax {
            if !self.step() {
                self.t = tmax;
            }
            r.push((self.t, self.species.clone()));
        }
        r
    }
    // pub fn advance_until(&mut self, tmax: f64) {
    //     while self.t < tmax {
    //         if !self.step() {
    //             self.t = tmax;
    //         }
    //     }
    // }
}

pub struct Gillespie<T: AsIndex> {
    species: Vec<isize>,
    t: f64,
    rates: Vec<Rate<T>>,
    reactions: Vec<Vec<isize>>
}

pub struct Rate<T: AsIndex> {
    rate: f64,
    species: Vec<SRate<T>>,
}

#[derive(Clone)]
pub enum SRate<T: AsIndex> {
    LMA(T), // S
    LMA2(T), // S
    MM(T, f64), // S, KC
    PosHill(T, f64, f64), // S, KC, n
    NegHill(T, f64, f64), // S, KC, n
}

impl<T: AsIndex+Clone> Rate<T> {
    pub fn new(rate: f64, species: &[SRate<T>]) -> Self {
        Rate {
            rate,
            species: species.to_vec()
        }
    }
    pub fn rate(&self, species: &[isize]) -> f64 {
        let mut r = self.rate;
        for factor in &self.species {
            match *factor {
                SRate::LMA(ref s) =>
                    r *= species[s.as_index()] as f64,
                SRate::LMA2(ref s) =>
                    r *= (species[s.as_index()]*species[s.as_index()]) as f64,
                SRate::MM(ref s, k) =>
                    r *= species[s.as_index()] as f64 /
                    (k + species[s.as_index()] as f64),
                SRate::PosHill(ref s, k, n) =>
                    r *= (species[s.as_index()] as f64).powf(n) /
                    (k.powf(n) + (species[s.as_index()] as f64).powf(n)),
                SRate::NegHill(ref s, k, n) =>
                    r *= (1. + (species[s.as_index()] as f64 / k).powf(n)).recip(),
            }
        }
        r
    }
}

#[cfg(test)]
mod tests {
    use gillespie::{Rate, SRate, Gillespie, AsIndex};
    #[test]
    fn sir() {
        index_enum! { enum SIR { S, I, R } }
        let mut sir = Gillespie::new(&[9999, 1, 0]);
        sir.add_reaction(Rate::new(0.1/10000., &[SRate::LMA(SIR::S), SRate::LMA(SIR::I)]), &[-1, 1, 0]);
        sir.add_reaction(Rate::new(0.01, &[SRate::LMA(SIR::I)]), &[0, -1, 1]);
        sir.advance_until(250.);
        assert_eq!(sir.get_species(&SIR::S) + sir.get_species(&SIR::I) + sir.get_species(&SIR::R), 10000);
    }
    #[test]
    fn dimers() {
        index_enum! { enum Dimers { G, M, P, D } }
        let mut dimers = Gillespie::new(&[1, 0, 0, 0]);
        dimers.add_reaction(Rate::new(25., &[SRate::LMA(Dimers::G)]), &[0, 1, 0, 0]);
        dimers.add_reaction(Rate::new(1000., &[SRate::LMA(Dimers::M)]), &[0, 0, 1, 0]);
        dimers.add_reaction(Rate::new(0.001, &[SRate::LMA2(Dimers::P)]), &[0, 0, -2, 1]);
        dimers.add_reaction(Rate::new(0.1, &[SRate::LMA(Dimers::M)]), &[0, -1, 0, 0]);
        dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::P)]), &[0, 0, -1, 0]);
        dimers.advance_until(1.);
        assert_eq!(dimers.get_species(&Dimers::G), 1);
        assert!(1000 < dimers.get_species(&Dimers::D));
        assert!(dimers.get_species(&Dimers::D) < 10000);
    }
}


#[cfg(test)]
mod benchmarks {
    use test::Bencher;
    use gillespie::{Rate, SRate, Gillespie, AsIndex};
    #[bench]
    fn sir(b: &mut Bencher) {
        index_enum! { enum SIR { S, I, R } }
        b.iter(|| {
            let mut sir = Gillespie::new(&[9999, 1, 0]);
            sir.add_reaction(Rate::new(0.1/10000., &[SRate::LMA(SIR::S), SRate::LMA(SIR::I)]), &[-1, 1, 0]);
            sir.add_reaction(Rate::new(0.05, &[SRate::LMA(SIR::R)]), &[0, -1, 1]);
            sir.advance_until(1000.);
        });
    }
    #[bench]
    fn dimers(b: &mut Bencher) {
        index_enum! { enum Dimers { G, M, P, D } }
        b.iter(|| {
            let mut dimers = Gillespie::new(&[1, 0, 0, 0]);
            dimers.add_reaction(Rate::new(25., &[SRate::LMA(Dimers::G)]), &[0, 1, 0, 0]);
            dimers.add_reaction(Rate::new(1000., &[SRate::LMA(Dimers::M)]), &[0, 0, 1, 0]);
            dimers.add_reaction(Rate::new(0.001, &[SRate::LMA2(Dimers::P)]), &[0, 0, -2, 1]);
            dimers.add_reaction(Rate::new(0.1, &[SRate::LMA(Dimers::M)]), &[0, -1, 0, 0]);
            dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::P)]), &[0, 0, -1, 0]);
            dimers.advance_until(1.);
        });
    }
    #[bench]
    fn dimers2(b: &mut Bencher) {
        index_enum! { enum Dimers { A, A_A, AA } }
        b.iter(|| {
            let mut dimers = Gillespie::new(&[100000, 0, 0]);
            dimers.add_reaction(Rate::new(1., &[SRate::LMA(Dimers::A)]), &[-1, 0, 0]);
            dimers.add_reaction(Rate::new(1./500., &[SRate::LMA2(Dimers::A)]), &[-2, 1, 0]);
            dimers.add_reaction(Rate::new(0.5, &[SRate::LMA(Dimers::A_A)]), &[2, -1, 0]);
            dimers.add_reaction(Rate::new(1./25., &[SRate::LMA(Dimers::A_A)]), &[0, -1, 1]);
            dimers.advance_until(30.);
        });
    }
}
