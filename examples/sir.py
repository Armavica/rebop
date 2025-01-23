"""The classical SIR model.

There are three categories of people:
* S: Susceptible
* I: Infected
* R: Recovered

    S + I -> 2 I
    I     ->   R
"""

import rebop

sir = rebop.Gillespie()
sir.add_reaction(1e-4, ["S", "I"], ["I", "I"])
sir.add_reaction(0.01, ["I"], ["R"])

print(sir)

ds = sir.run({"S": 999, "I": 1}, tmax=250, nb_steps=250, rng=0)

print(ds)
