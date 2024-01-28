import rebop

alphaA = 50
alphapA = 500
alphaR = 0.01
alphapR = 50
betaA = 50
betaR = 5
deltaMA = 10
deltaMR = 0.5
deltaA = 1
deltaR = 0.2
gammaA = 1
gammaR = 1
gammaC = 2
thetaA = 50
thetaR = 100

vilar = rebop.Gillespie()
vilar.add_reaction(gammaA, ['Da', 'A'], ['Dpa'])
vilar.add_reaction(gammaR, ['Dr', 'A'], ['Dpr'])
vilar.add_reaction(thetaA, ['Dpa'], ['Da', 'A'])
vilar.add_reaction(thetaR, ['Dpr'], ['Dr', 'A'])
vilar.add_reaction(alphaA, ['Da'], ['Da', 'Ma'])
vilar.add_reaction(alphaR, ['Dr'], ['Dr', 'Mr'])
vilar.add_reaction(alphapA, ['Dpa'], ['Dpa', 'Ma'])
vilar.add_reaction(alphapR, ['Dpr'], ['Dpr', 'Mr'])
vilar.add_reaction(betaA, ['Ma'], ['Ma', 'A'])
vilar.add_reaction(betaR, ['Mr'], ['Mr', 'R'])
vilar.add_reaction(gammaC, ['A', 'R'], ['C'])
vilar.add_reaction(deltaA, ['C'], ['R'])
vilar.add_reaction(deltaMA, ['Ma'], [])
vilar.add_reaction(deltaMR, ['Mr'], [])
vilar.add_reaction(deltaA, ['A'], [])
vilar.add_reaction(deltaR, ['R'], [])

print(vilar)

for _ in range(100):
    ds = vilar.run({'Da': 1, 'Dr': 1}, tmax=200, nb_steps=200)
