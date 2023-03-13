# Vilar oscillator benchmark

We simulated the Vilar oscillator model with several programs.

![Vilar oscillator benchmark](vilar.png)

## rebop (Rust)

`hyperfine 'touch vilar.rs && cargo run --release --bin vilar'`

## rebopy (Python bindings to rebop)

`cargo build --release && ln -s ../../target/release/librebop.so rebop.so`

`time python3 vilar.py`

## GillesPy2 (C++ & Python)

`time python3 vilar_gillespy2.py`

## BioNetGen (C++ & Perl)

`time BNG2.pl vilar.bngl`

## DiffEqJump.jl (Julia)

`time julia vilar.jl`
