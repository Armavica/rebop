//! Bebop is a fast stochastic simulator for chemical reaction networks.
//!
//! It offers two independent APIs to describe and simulate reaction
//! networks: a macro-based implemented by [`define_system`], and a
//! function-based implemented by the module [`gillespie`].
//!

pub mod gillespie;
mod gillespie_macro;
