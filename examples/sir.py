import rebop

sir = rebop.Gillespie()
sir.add_reaction(1e-4, ["S", "I"], ["I", "I"])
sir.add_reaction(0.01, ["I"], ["R"])
print(sir)

df = sir.run({"S": 999, "I": 1}, tmax=250, nb_steps=250)

print(df)
