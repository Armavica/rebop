## rebop

rebop is a stochastic simulator for chemical reaction networks, written
in Rust and designed to be as fast and convenient as possible.

### Example

Consider the following dimerization system:

* a mRNA can be transcribed from a gene;
* a protein can be produced from a mRNA;
* two proteins can dimerize;
* mRNA and proteins are subject to constant decay;

To simulate this system, we can write the following:

``` rust
use rebop::define_system;

define_system! {
    r_tx r_tl r_dim r_decay_mrna r_decay_prot;
    Dimers { gene, mRNA, protein, dimer }
    transcription   : gene              => gene, mRNA       @ r_tx
    translation     : mRNA              => mRNA, protein    @ r_tl
    dimerization    : protein, protein  => dimer            @ r_dim
    decay_mRNA      : mRNA              =>                  @ r_decay_mrna
    decay_protein   : protein           =>                  @ r_decay_prot
}

fn main() {
    let mut dimers = Dimers::new();
    dimers.r_tx = 25.;
    dimers.r_tl = 1000.;
    dimers.r_dim = 0.001;
    dimers.r_decay_mrna = 0.1;
    dimers.r_decay_prot = 1.;
    dimers.gene = 1;
    dimers.advance_until(1.);
    println!("t = {}: dimer = {}", dimers.t, dimers.dimer);
}
```


