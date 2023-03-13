from cayenne.simulation import Simulation
# import matplotlib.pyplot as plt

model_str = """
const compartment comp1;
comp1 = 1.0; # volume of compartment

activation_a        : Da + A    => Dpa      ; gammaA;
activation_r        : Dr + A    => Dpr      ; gammaR;
deactivation_a      : Dpa       => Da + A   ; thetaA;
deactivation_r      : Dpr       => Dr + A   ; thetaR;
transcription_a     : Da        => Da + Ma  ; alphaA;
transcription_r     : Dr        => Dr + Mr  ; alphaR;
transcription_p_a   : Dpa       => Dpa + Ma ; alphapA;
transcription_p_r   : Dpr       => Dpr + Mr ; alphapR;
translation_a       : Ma        => Ma + A   ; betaA;
translation_r       : Mr        => Mr + R   ; betaR;
complexation        : A + R     => C        ; gammaC;
decomplexation      : C         => R        ; gammaA;
decay_mRNA_a        : Ma        =>          ; deltaMA;
decay_mRNA_r        : Mr        =>          ; deltaMR;
decay_prot_a        : A         =>          ; deltaA;
decay_prot_r        : R         =>          ; deltaR;

alphaA = 50;
alphapA = 500;
alphaR = 0.01;
alphapR = 50;
betaA = 50;
betaR = 5;
deltaMA = 10;
deltaMR = 0.5;
deltaA = 1;
deltaR = 0.2;
gammaA = 1;
gammaR = 2;
gammaC = 2;
thetaA = 50;
thetaR = 100;
chem_flag = false;

Da = 1;
Dr = 1;
Dpa = 0;
Dpr = 0;
Ma = 0;
Mr = 0;
A = 0;
R = 0;
C = 0;
"""
sim = Simulation.load_model(model_str, "ModelString")
# Run the simulation
sim.simulate(max_t=200, max_iter=700000, algorithm='direct')
print(sim.results)
# sim.plot()
# plt.show()
