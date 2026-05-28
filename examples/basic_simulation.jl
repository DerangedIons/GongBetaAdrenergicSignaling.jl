"""
Basic simulation example for the Gong beta-adrenergic signaling model.

Demonstrates the plain-ODE path: build the parameter vector, solve `rhs_signaling!`
with OrdinaryDiffEq.jl, and compute phosphorylation fractions from the solution.
"""

using Pkg
Pkg.activate(@__DIR__)
using GongBetaAdrenergicSignaling
using OrdinaryDiffEq

N_beats = 200
BCL = 1000

# 1 μM isoproterenol stimulation
p = compute_parameters(1.0)
u0 = default_initial_state()
prob = ODEProblem(rhs_signaling!, u0, (0.0, N_beats * BCL), p)
sol = solve(prob, Tsit5())

# Phosphorylation fractions at the start and end of the simulation
labels = ("fICaL", "fIKs", "fPLB", "fTnI", "fINa", "fINaK", "fRyR", "fIKur")
f_start = zeros(8)
f_end = zeros(8)
effective_fractions!(f_start, sol.u[1], p)
effective_fractions!(f_end, sol.u[end], p)

println("\n=== Phosphorylation fractions (t=0 → t=$(sol.t[end]) ms) ===")
for (i, label) in enumerate(labels)
    println("  $label:  ", f_start[i], " → ", f_end[i])
end
