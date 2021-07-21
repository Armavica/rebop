# bebop

Two goals of this project are efficiency and practicality.  The
following macro defines the corresponding reaction network:

``` Rust
define_system! {
    Dimers { gene, mRNA, protein, dimer }
    r_transcription : gene             => gene, mRNA    @ 25.
    r_translation   : mRNA             => mRNA, protein @ 1000.
    r_dimerization  : protein, protein => dimer         @ 0.001
    r_decay_mRNA    : mRNA             => nil           @ 0.1
    r_decay_protein : protein          => nil           @ 1.
}
```

To simulate the system, instantiate a new problem, set the initial
values, and simulate:

``` Rust
let mut problem = Dimers::new();
problem.gene = 1;
let trace = problem.advance_until(1.);
println!("{}: dimer = {}", problem.t, problem.dimer);
```

Or for the classic SIR example:

``` Rust
define_system! {
    SIR { S, I, R }
    r_infection: S, I => I, I @ 0.1/1000.
    r_healing  : I    => R    @ 0.01
}

fn main() {
    let mut problem = SIR::new();
    problem.S = 999;
    problem.I = 1;
    let trace = problem.advance_until(250.);
    println!("time,S,I,R");
    for state in &trace {
        println!("{},{},{},{}", state.t, state.S, state.I, state.R);
    }
}
```

which can produce an output such as
![SIR](https://github.com/Armavica/bebop/blob/master/sir.png)

## Not (yet) features

* propensities != reaction rates
* Next reaction method (Gibson - Brooke)
* Tau-leaping
* Adaptive tau-leaping
* Hybrid stoch / diff
* Space (localizations)
* Diffusion
* Volume change
* Michaelis-Menten
* Time-varying inputs
* Rule modeling
* SBML
* Parameter estimation
* Local sensitivity analysis
* Parallelization

### Benchmarks

* Dimers
* SIR
* Toggle switch
* STOCKS
    * LacZ (example 1)
    * LacY/LacZ (example 2)
* StochSim
    * Lotka Volterra
    * Michaelis-Menten
    * Network
* SimBiology
    * G Protein
* Cellware
    * Brusselator / Oregonator
* Dizzy
    * GAL
    * Repressilator

### Other software

* [GillesPy2](https://github.com/StochSS/GillesPy2)
* [STEPS](https://github.com/CNS-OIST/STEPS)
* [SimBiology](https://fr.mathworks.com/help/simbio/)
* [Copasi](http://copasi.org/)
* [BioNetGen](http://bionetgen.org/)
* [VCell](http://vcell.org/)
* [Smoldyn](http://www.smoldyn.org/)
* [KaSim](https://kappalanguage.org/)

#### Seem unmaintained

* [Dizzy](http://magnet.systemsbiology.net/software/Dizzy/)
* [Cellware](http://www.bii.a-star.edu.sg/achievements/applications/cellware/)
* [STOCKS](https://doi.org/10.1093/bioinformatics/18.3.470)
* [StochSim](http://lenoverelab.org/perso/lenov/stochsim.html)
* [Systems biology toolbox](http://www.sbtoolbox.org/)
* [StochKit](https://github.com/StochSS/StochKit)
* [SmartCell](http://software.crg.es/smartcell/)
* [NFsim](http://michaelsneddon.net/nfsim/)

