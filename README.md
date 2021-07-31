# bebop

Bebop is a fast stochastic simulator for well-mixed chemical reaction
networks.

Two goals of this project are efficiency and convenience.
The following macro defines a reaction network in a natural manner:

``` rust
define_system! {
    r_tx r_tl r_dim r_decay_mrna r_decay_prot;
    Dimers { gene, mRNA, protein, dimer }
    transcription : gene             => gene, mRNA    @ r_tx
    translation   : mRNA             => mRNA, protein @ r_tl
    dimerization  : protein, protein => dimer         @ r_dim
    decay_mRNA    : mRNA             =>               @ r_decay_mrna
    decay_protein : protein          =>               @ r_decay_prot
}
```

To simulate the system: instantiate a new problem, set the initial
values, the parameters, and launch the simulation.

``` rust
let mut problem = Dimers::with_parameters(25., 1000., 0.001, 0.1, 1.);
problem.gene = 1;
problem.advance_until(1.);
println!("{}: dimer = {}", problem.t, problem.dimer);
```

Or for the classic SIR example:

``` rust
define_system! {
    r_inf r_heal;
    SIR { S, I, R }
    infection: S, I => I, I @ r_inf
    healing  : I    => R    @ r_heal
}

fn main() {
    let mut problem = SIR::new();
    problem.r_inf = 0.1 / 1000.;
    problem.r_heal = 0.01;
    problem.S = 999;
    problem.I = 1;
    println!("time,S,I,R");
    for t in 0..250 {
        problem.advance_until(t as f64);
        println!("{},{},{},{}", problem.t, problem.S, problem.I, problem.R);
    }
}
```

which can produce an output such as
![SIR](https://github.com/Armavica/bebop/blob/master/sir.png)

## Performance

On typical example networks, bebop outperformed all other software.

*Disclaimer*: Most of these softwares contain much more features than
bebop (e.g. spatial models, custom reaction rates, etc.).  Some of these
features might require them to make compromises on speed.  Moreover,
some can be conveniently used through wrappers (for example when the
simulation code is written in C++ but the model is expressible in
Python).  These wrappers can also add a significant overhead.

To benchmark these softwares in the fairest possible conditions, we
considered for everyone a typical situation where the model was just
modified and we want to simulate it `N` times.  So the (re-)compilation
time is included in this benchmark.

Example for the Vilar oscillator (*Mechanisms of noise-resistance in
genetic oscillators*, Vilar et al., PNAS 2002).  Here, we simulate this
model from `t=0` to `t=200`, saving the state at time intervals of `1`.

![Vilar oscillator performance](https://github.com/Armavica/bebop/blob/master/vilar.png)

Bebop is the fastest, both per simulation, and with compilation time
included.

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

