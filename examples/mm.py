"""Demonstration of a non LMA rate and of rate parameters.

The reaction here is A → B with the enzymatic rate V * A / (Km + A).

V and Km are provided as parameters in the `run` method.
"""

import rebop

mm = rebop.Gillespie()
mm.add_reaction("V * A / (Km + A)", ["A"], ["P"])

print(mm)

params = {"V": 1, "Km": 20}

ds = mm.run({"A": 100}, tmax=250, nb_steps=100, params=params, rng=0)

print(ds)
