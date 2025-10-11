use pyo3::exceptions::{PyUserWarning, PyValueError};
use pyo3::prelude::*;

use std::collections::HashMap;
use std::fmt::{self, Display, Formatter};

pub use rand;
pub use rand_distr;

use crate::expr::{PExpr, RateParseError};
use crate::gillespie;

/// Reaction system composed of species and reactions.
#[pyclass]
struct Gillespie {
    species: HashMap<String, usize>,
    init: HashMap<String, usize>,
    reactions: Vec<(PRate, Vec<String>, Vec<String>)>,
}

#[derive(Clone, Debug, FromPyObject)]
enum PyRate {
    Lma(f64),
    Expr(String),
}

#[derive(Clone, Debug)]
enum PRate {
    Lma(f64, Vec<String>),
    Expr(PExpr),
}

impl Display for PRate {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match self {
            PRate::Lma(rate, _) => write!(f, "LMA({rate})"),
            PRate::Expr(expr) => write!(f, "{expr}"),
        }
    }
}

impl PRate {
    fn new(rate: PyRate, reactants: &[String]) -> Result<Self, RateParseError> {
        match rate {
            PyRate::Lma(c) => Ok(PRate::Lma(c, reactants.to_vec())),
            PyRate::Expr(rate) => rate.parse().map(PRate::Expr),
        }
    }
    fn to_gillespie_rate(
        &self,
        species: &HashMap<String, usize>,
        params: &HashMap<String, f64>,
    ) -> Result<gillespie::Rate, String> {
        let rate = match self {
            PRate::Lma(rate, reactants) => {
                let mut rate_reactants = vec![0; species.len()];
                for reactant in reactants {
                    rate_reactants[species[reactant]] += 1;
                }
                gillespie::Rate::lma(*rate, rate_reactants)
            }
            PRate::Expr(pexpr) => gillespie::Rate::expr(pexpr.to_expr(species, params)?),
        };
        Ok(rate)
    }
}

impl Gillespie {
    fn add_species_from_expr(&mut self, expr: &PExpr) {
        match expr {
            PExpr::Constant(_) => {}
            PExpr::Variable(s) => self.add_species(s),
            PExpr::Add(a, b)
            | PExpr::Sub(a, b)
            | PExpr::Mul(a, b)
            | PExpr::Div(a, b)
            | PExpr::Pow(a, b) => {
                self.add_species_from_expr(a);
                self.add_species_from_expr(b);
            }
            PExpr::Exp(a) => {
                self.add_species_from_expr(a);
            }
        }
    }
}

#[pymethods]
impl Gillespie {
    #[new]
    fn new() -> Self {
        Gillespie {
            species: HashMap::new(),
            init: HashMap::new(),
            reactions: Vec::new(),
        }
    }
    /// Register a species to the model.
    fn add_species(&mut self, species: &str) {
        if !self.species.contains_key(species) {
            self.species.insert(species.to_string(), self.species.len());
        }
    }
    /// Number of species currently in the system.
    fn nb_species(&self) -> PyResult<usize> {
        Ok(self.species.len())
    }
    /// Add a reaction to the system.
    #[pyo3(signature = (rate, reactants, products, reverse_rate=None))]
    fn add_reaction(
        &mut self,
        rate: PyRate,
        reactants: Vec<String>,
        products: Vec<String>,
        reverse_rate: Option<PyRate>,
    ) -> PyResult<()> {
        // Convert PyRate into PRate (and possibly fail with rate parse error)
        let rate = match PRate::new(rate, &reactants) {
            Ok(rate) => rate,
            Err(_) => return Err(PyValueError::new_err("Rate expression not understood")),
        };
        // Insert unknown reactants in known species
        for reactant in &reactants {
            self.add_species(reactant);
        }
        // Insert unknown products in known species
        for product in &products {
            self.add_species(product);
        }
        self.reactions
            .push((rate, reactants.clone(), products.clone()));
        if let Some(rrate) = reverse_rate {
            self.add_reaction(rrate, products, reactants, None)?;
        }
        Ok(())
    }
    /// Number of reactions currently in the system.
    fn nb_reactions(&self) -> PyResult<usize> {
        Ok(self.reactions.len())
    }
    /// Set the initial count of species
    fn set_init(&mut self, init: HashMap<String, usize>) -> PyResult<()> {
        let mut warning = false;
        for species in init.keys() {
            if !self.species.contains_key(species) {
                warning = true;
                self.add_species(species);
            }
        }
        self.init = init;
        if warning {
            let message = "Some species are not involved in any reactions. You should probably instead use parameters.";
            Err(PyUserWarning::new_err(message))
        } else {
            Ok(())
        }
    }
    /// Run the system until `tmax` with `nb_steps` steps.
    ///
    /// The initial configuration is specified in the dictionary `init`.
    /// Returns `times, vars` where `times` is an array of `nb_steps + 1` uniformly spaced time
    /// points between `0` and `tmax`, and `vars` is a dictionary of species name to array of
    /// values at the given time points.  One can specify a random `seed` for reproducibility.
    /// If `nb_steps` is `0`, then returns all reactions, ending with the first that happens at
    /// or after `tmax`.
    #[pyo3(signature = (tmax, nb_steps, params, seed=None, sparse=None, var_names=None))]
    fn run(
        &self,
        tmax: f64,
        nb_steps: usize,
        params: HashMap<String, f64>,
        seed: Option<u64>,
        sparse: Option<bool>,
        var_names: Option<Vec<String>>,
    ) -> PyResult<(Vec<f64>, HashMap<String, Vec<isize>>)> {
        for p in params.keys() {
            if self.species.contains_key(p) {
                return Err(PyValueError::new_err(format!(
                    "Species {p} cannot also be a parameter."
                )));
            }
        }
        let sparse = sparse.unwrap_or(false);
        let mut x0 = vec![0; self.species.len()];
        for (name, &value) in &self.init {
            if let Some(&id) = self.species.get(name) {
                x0[id] = value as isize;
            }
        }
        let mut g = match seed {
            Some(seed) => gillespie::Gillespie::new_with_seed(x0, sparse, seed),
            None => gillespie::Gillespie::new(x0, sparse),
        };
        let save_indices: Vec<_> = match &var_names {
            Some(x) => x
                .iter()
                .map(|key| *self.species.get(key).unwrap())
                .collect(),
            None => (0..self.species.len()).collect(),
        };

        for (rate, reactants, products) in self.reactions.iter() {
            let mut vreactants = vec![0; self.species.len()];
            for reactant in reactants {
                vreactants[self.species[reactant]] += 1;
            }
            let mut actions = vec![0; self.species.len()];
            for reactant in reactants {
                actions[self.species[reactant]] -= 1;
            }
            for product in products {
                actions[self.species[product]] += 1;
            }
            match rate.to_gillespie_rate(&self.species, &params) {
                Ok(rate) => g.add_reaction(rate, actions),
                Err(e) => return Err(PyValueError::new_err(e)),
            }
        }
        let mut times = Vec::new();
        // species.shape = (species, nb_steps)
        let mut species = vec![Vec::new(); save_indices.len()];
        if nb_steps > 0 {
            for i in 0..=nb_steps {
                let t = tmax * i as f64 / nb_steps as f64;
                times.push(t);
                g.advance_until(t);
                for (i, s) in save_indices.iter().enumerate() {
                    species[i].push(g.get_species(*s));
                }
            }
        } else {
            // nb_steps = 0: we return every step
            let mut rates = vec![f64::NAN; g.nb_reactions()];
            times.push(g.get_time());
            for (i, s) in save_indices.iter().enumerate() {
                species[i].push(g.get_species(*s));
            }
            while g.get_time() < tmax {
                g._advance_one_reaction(&mut rates);
                times.push(g.get_time());
                for (i, s) in save_indices.iter().enumerate() {
                    species[i].push(g.get_species(*s));
                }
            }
        }
        let mut result = HashMap::new();
        match var_names {
            Some(x) => {
                for (id, name) in x.iter().enumerate() {
                    result.insert(name.clone(), species[id].clone());
                }
            }
            None => {
                for (name, &id) in &self.species {
                    result.insert(name.clone(), species[id].clone());
                }
            }
        }
        Ok((times, result))
    }
    fn __str__(&self) -> PyResult<String> {
        let mut s = format!(
            "{} species and {} reactions\n",
            self.species.len(),
            self.reactions.len()
        );
        for (rate, reactants, products) in &self.reactions {
            s.push_str(&reactants.join(" + "));
            s.push_str(" --> ");
            s.push_str(&products.join(" + "));
            s.push_str(&format!(" @ {rate}\n"));
        }
        Ok(s)
    }
}

#[pymodule]
fn _lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<Gillespie>()?;
    Ok(())
}
