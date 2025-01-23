//! rebop is a fast stochastic simulator for well-mixed chemical
//! reaction networks.
//!
//! Performance and ergonomics are taken very seriously.  For this reason,
//! two independent APIs are provided to describe and simulate reaction
//! networks:
//!
//! * a macro-based DSL implemented by [`define_system`], usually the
//! most efficient, but that requires to compile a rust program;
//! * a function-based API implemented by the module [`gillespie`], also
//! available through Python bindings.  This one does not require a rust
//! compilation and allows the system to be defined at run time.  It is
//! typically 2 or 3 times slower than the macro DSL, but still faster
//! than all other software tried.
//!
//! # The macro DSL
//!
//! It currently only supports reaction rates defined by the law of mass
//! action.  The following macro defines a dimerization reaction network
//! naturally:
//!
//! ```rust
//! use rebop::define_system;
//! define_system! {
//!     r_tx r_tl r_dim r_decay_mRNA r_decay_prot;
//!     Dimers { gene, mRNA, protein, dimer }
//!     transcription   : gene      => gene + mRNA      @ r_tx
//!     translation     : mRNA      => mRNA + protein   @ r_tl
//!     dimerization    : 2 protein => dimer            @ r_dim
//!     decay_mRNA      : mRNA      =>                  @ r_decay_mRNA
//!     decay_protein   : protein   =>                  @ r_decay_prot
//! }
//! ```
//!
//! To simulate the system, put this definition in a rust code file and
//! instantiate the problem, set the parameters, the initial values, and
//! launch the simulation:
//!
//! ```rust
//! # use rebop::define_system;
//! # define_system! {
//! #     r_tx r_tl r_dim r_decay_mRNA r_decay_prot;
//! #     Dimers { gene, mRNA, protein, dimer }
//! #     transcription   : gene      => gene + mRNA      @ r_tx
//! #     translation     : mRNA      => mRNA + protein   @ r_tl
//! #     dimerization    : 2 protein => dimer            @ r_dim
//! #     decay_mRNA      : mRNA      =>                  @ r_decay_mRNA
//! #     decay_protein   : protein   =>                  @ r_decay_prot
//! # }
//! let mut problem = Dimers::new();
//! problem.r_tx = 25.0;
//! problem.r_tl = 1000.0;
//! problem.r_dim = 0.001;
//! problem.r_decay_mRNA = 0.1;
//! problem.r_decay_prot = 1.0;
//! problem.gene = 1;
//! problem.advance_until(1.0);
//! println!("t = {}: dimer = {}", problem.t, problem.dimer);
//! ```
//!
//! Or for the classic SIR example:
//!
//! ```rust
//! use rebop::define_system;
//!
//! define_system! {
//!     r_inf r_heal;
//!     SIR { S, I, R }
//!     infection   : S + I => 2 I  @ r_inf
//!     healing     : I     => R    @ r_heal
//! }
//!
//! fn main() {
//!     let mut problem = SIR::new();
//!     problem.r_inf = 1e-4;
//!     problem.r_heal = 0.01;
//!     problem.S = 999;
//!     problem.I = 1;
//!     println!("time,S,I,R");
//!     for t in 0..250 {
//!         problem.advance_until(t as f64);
//!         println!("{},{},{},{}", problem.t, problem.S, problem.I, problem.R);
//!     }
//! }
//! ```
//!
//! which can produce an output similar to this one:
//!
//! ![Typical SIR output](https://github.com/Armavica/rebop/blob/main/sir.png?raw=true)
//!
//! # Python bindings
//!
//! This API shines through the Python bindings which allow one to
//! define a model easily:
//!
//! ```python
//! import rebop
//!
//! sir = rebop.Gillespie()
//! sir.add_reaction(1e-4, ['S', 'I'], ['I', 'I'])
//! sir.add_reaction(0.01, ['I'], ['R'])
//! print(sir)
//!
//! ds = sir.run({'S': 999, 'I': 1}, tmax=250, nb_steps=250)
//! ```
//!
//! You can test this code by installing `rebop` from PyPI with
//! `pip install rebop`. To build the Python bindings from source,
//! the simplest is to clone this git repository and use `maturin
//! develop`.
//!
//! # The traditional API
//!
//! The function-based API underlying the Python package is also available
//! from Rust, if you want to be able to define models at run time (instead
//! of at compilation time with the macro DSL demonstrated above).
//! The SIR model is defined as:
//!
//! ```rust
//! use rebop::gillespie::{Gillespie, Rate};
//!
//! let mut sir = Gillespie::new([999, 1, 0], false);
//! //                           [  S, I, R]
//! // S + I => 2 I with rate 1e-4
//! sir.add_reaction(Rate::lma(1e-4, [1, 1, 0]), [-1, 1, 0]);
//! // I => R with rate 0.01
//! sir.add_reaction(Rate::lma(0.01, [0, 1, 0]), [0, -1, 1]);
//!
//! println!("time,S,I,R");
//! for t in 0..250 {
//!     sir.advance_until(t as f64);
//!     println!("{},{},{},{}", sir.get_time(), sir.get_species(0), sir.get_species(1), sir.get_species(2));
//! }
//! ```
//!
//! # Performance
//!
//! Performance is taken very seriously, and as a result, rebop
//! outperforms every other package and programming language that we
//! tried.
//!
//! *Disclaimer*: Most of this software currently contains much more
//! features than rebop (e.g. spatial models, custom reaction rates,
//! etc.).  Some of these features might have required them to make
//! compromises on speed.  Moreover, as much as we tried to keep the
//! comparison fair, some return too much or too little data, or write
//! them on disk.  The baseline that we tried to approach for all these
//! programs is the following: *the model was just modified, we want
//! to simulate it `N` times and print regularly spaced measurement
//! points*.  This means that we always include initialization or
//! (re-)compilation time if applicable.  We think that it is the most
//! typical use-case of a researcher who works on the model.  This
//! benchmark methods allows to record both the initialization time
//! (y-intercept) and the simulation time per simulation (slope).
//!
//! Many small benchmarks on toy examples are tracked to guide the
//! development.  To compare the performance with other software,
//! we used a real-world model of low-medium size (9 species and 16
//! reactions): the Vilar oscillator (*Mechanisms of noise-resistance
//! in genetic oscillators*, Vilar et al., PNAS 2002).  Here, we
//! simulate this model from `t=0` to `t=200`, reporting the state at
//! time intervals of `1` time unit.
//!
//! ![Vilar oscillator benchmark](https://github.com/Armavica/rebop/blob/main/benches/vilar/vilar.png?raw=true)
//!
//! We can see that rebop's macro DSL is the fastest of all, both in
//! time per simulation, and with compilation time included.  The second
//! fastest is rebop's traditional API invoked by convenience through
//! the Python bindings.
//!
//! # Features to come
//!
//! * compartment volumes
//! * arbitrary reaction rates
//! * other SSA algorithms
//! * tau-leaping
//! * adaptive tau-leaping
//! * hybrid models (continuous and discrete)
//! * SBML
//! * CLI interface
//! * parameter estimation
//! * local sensitivity analysis
//! * parallelization
//!
//! # Features probably not to come
//!
//! * events
//! * space (reaction-diffusion systems)
//! * rule modelling
//!
//! # Benchmark ideas
//!
//! * DSMTS
//! * purely decoupled exponentials
//! * ring
//! * Toggle switch
//! * LacZ, LacY/LacZ (from STOCKS)
//! * Lotka Volterra, Michaelis--Menten, Network (from StochSim)
//! * G protein (from SimBiology)
//! * Brusselator / Oregonator (from Cellware)
//! * GAL, repressilator (from Dizzy)
//!
//! # Similar software
//!
//! ## Maintained
//!
//! * [GillesPy2](https://github.com/StochSS/GillesPy2)
//! * [STEPS](https://github.com/CNS-OIST/STEPS)
//! * [SimBiology](https://fr.mathworks.com/help/simbio/)
//! * [Copasi](http://copasi.org/)
//! * [BioNetGen](http://bionetgen.org/)
//! * [VCell](http://vcell.org/)
//! * [Smoldyn](http://www.smoldyn.org/)
//! * [KaSim](https://kappalanguage.org/)
//! * [StochPy](https://github.com/SystemsBioinformatics/stochpy)
//! * [BioSimulator.jl](https://github.com/alanderos91/BioSimulator.jl)
//! * [DiffEqJump.jl](https://github.com/SciML/DiffEqJump.jl)
//! * [Gillespie.jl](https://github.com/sdwfrost/Gillespie.jl)
//! * [GillespieSSA2](https://github.com/rcannood/GillespieSSA2)
//! * [Cayenne](https://github.com/quantumbrake/cayenne)
//!
//! ## Seem unmaintained
//!
//! * [Dizzy](http://magnet.systemsbiology.net/software/Dizzy/)
//! * [Cellware](http://www.bii.a-star.edu.sg/achievements/applications/cellware/)
//! * [STOCKS](https://doi.org/10.1093/bioinformatics/18.3.470)
//! * [StochSim](http://lenoverelab.org/perso/lenov/stochsim.html)
//! * [Systems biology toolbox](http://www.sbtoolbox.org/)
//! * [StochKit](https://github.com/StochSS/StochKit) (successor: GillesPy2)
//! * [SmartCell](http://software.crg.es/smartcell/)
//! * [NFsim](http://michaelsneddon.net/nfsim/)

use pyo3::exceptions::{PyUserWarning, PyValueError};
use pyo3::prelude::*;
use std::collections::HashMap;
use std::fmt::{self, Display, Formatter};

pub use rand;
pub use rand_distr;

mod expr;
pub mod gillespie;
mod gillespie_macro;

use expr::PExpr;

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
    fn new(rate: PyRate, reactants: &[String]) -> Result<Self, expr::RateParseError> {
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
            let message =
                "Some species are not involved in any reactions. You should probably instead use parameters.";
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
