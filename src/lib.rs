//! rebop is a fast stochastic simulator for chemical reaction networks.
//!
//! It offers two independent APIs to describe and simulate reaction
//! networks: a macro-based implemented by [`define_system`], and a
//! function-based implemented by the module [`gillespie`].
//!

use pyo3::prelude::*;
use pyo3::class::basic::PyObjectProtocol;
use std::collections::HashMap;

pub mod gillespie;
mod gillespie_macro;

#[pyclass]
struct Gillespie {
    species: HashMap<String, usize>,
    reactions: Vec<(f64, Vec<String>, Vec<String>)>,
}

#[pymethods]
impl Gillespie {
    #[new]
    fn new() -> Self {
        Gillespie {
            species: HashMap::new(),
            reactions: Vec::new(),
        }
    }
    fn nb_species(&self) -> PyResult<usize> {
        Ok(self.species.len())
    }
    fn add_reaction(&mut self, rate: f64, reactants: Vec<String>, products: Vec<String>) -> PyResult<()> {
        for reactant in &reactants {
            if !self.species.contains_key(reactant) {
                self.species.insert(reactant.clone(), self.species.len());
            }
        }
        for product in &products {
            if !self.species.contains_key(product) {
                self.species.insert(product.clone(), self.species.len());
            }
        }
        self.reactions.push((rate, reactants, products));
        Ok(())
    }
    fn nb_reactions(&self) -> PyResult<usize> {
        Ok(self.reactions.len())
    }
    fn run(&self, init: HashMap<String, usize>, tmax: f64, nb_steps: usize) -> PyResult<(Vec<f64>, HashMap<String, Vec<isize>>)> {
        let mut x0 = vec![0; self.species.len()];
        for (name, &value) in &init {
            if let Some(&id) = self.species.get(name) {
                x0[id] = value as isize;
            }
        }
        let mut g = gillespie::Gillespie::<usize>::new(x0);
        for (rate, reactants, products) in self.reactions.iter() {
            let rate = gillespie::Rate::new(
                *rate,
                &reactants.iter().map(|r| gillespie::SRate::LMA(self.species[r])).collect::<Vec<_>>());
            let mut actions = vec![0; self.species.len()];
            for reactant in reactants {
                actions[self.species[reactant]] -= 1;
            }
            for product in products {
                actions[self.species[product]] += 1;
            }
            g.add_reaction(rate, actions);
        }
        let mut times = Vec::new();
        let mut species = vec![Vec::new(); self.species.len()];
        for i in 0..=nb_steps {
            let t = tmax * i as f64 / nb_steps as f64;
            times.push(t);
            g.advance_until(t);
            for s in 0..self.species.len() {
                species[s].push(g.get_species(&s));
            }
        }
        let mut result = HashMap::new();
        for (name, &id) in &self.species {
            result.insert(name.clone(), species[id].clone());
        }
        Ok((times, result))
    }
}

#[pyproto]
impl<'p> PyObjectProtocol<'p> for Gillespie {
    fn __str__(&self) -> PyResult<String>{
        let mut s = String::new();
        s.push_str(&format!("{} species and {} reactions\n", self.species.len(), self.reactions.len()));
        for (rate, reactants, products) in &self.reactions {
            s.push_str(&reactants.join(" + "));
            s.push_str(" --> ");
            s.push_str(&products.join(" + "));
            s.push_str(&format!(" @ {}\n", rate));
        }
        Ok(s)
    }
}

#[pymodule]
fn rebop(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Gillespie>()?;
    Ok(())
}

