using Catalyst
using JumpProcesses

rn = @reaction_network begin
    @species Da(t) Dr(t) Dpa(t) Dpr(t) Ma(t) Mr(t) A(t) R(t) C(t)
    @parameters αA αpA αR αpR βA βR δMA δMR δA δR γA γR γC θA θR
    (γA, θA),   Da + A ↔ Dpa
    (γR, θR),   Dr + A ↔ Dpr
    αA,         Da → Da + Ma
    αR,         Dr → Dr + Mr
    αpA,        Dpa → Dpa + Ma
    αpR,        Dpr → Dpr + Mr
    βA,         Ma → Ma + A
    βR,         Mr → Mr + R
    γC,         A + R → C
    δA,         C → R
    δMA,        Ma → ∅
    δMR,        Mr → ∅
    δA,         A → ∅
    δR,         R → ∅
end

p = [50., 500., 0.01, 50., 50., 5., 10., 0.5, 1., 0.2, 1., 1., 2., 50., 100.]
u0 = [1, 1, 0, 0, 0, 0, 0, 0, 0]
tspan = (0., 200.)
prob_discrete = DiscreteProblem(rn, u0, tspan, p)
prob_jump = JumpProblem(rn, prob_discrete, Direct(), save_positions=(false, false))
for i = 1:100
    sol = solve(prob_jump, SSAStepper(), saveat=1.)
end
