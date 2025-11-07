"""
Basic simulation example for the Gong Beta Adrenergic Signaling model.

This example demonstrates:
1. Baseline simulation (no isoproterenol)
2. Printing phosphorylation observables before and after simulation
"""

using Pkg; Pkg.activate(@__DIR__)
using GongBetaAdrenergicSignaling
using OrdinaryDiffEq

N_beats = 200
BCL = 1000

# Baseline simulation
@mtkbuild sys = GongBetaAdrenergic(iso_conc=1.0)
prob = ODEProblem(sys, [], (0.0, N_beats * BCL))
sol = solve(prob, Tsit5())

# Print phosphorylation observables
println("\n=== Phosphorylation Observables (t=0 → t=$(sol.t[end]) ms) ===")
println("  fICaL_PKA:  ", sol[sys.fICaL_PKA][1], " → ", sol[sys.fICaL_PKA][end])
println("  fIKs_PKA:   ", sol[sys.fIKs_PKA][1], " → ", sol[sys.fIKs_PKA][end])
println("  fPLB_PKA:   ", sol[sys.fPLB_PKA][1], " → ", sol[sys.fPLB_PKA][end])
println("  fTnI_PKA:   ", sol[sys.fTnI_PKA][1], " → ", sol[sys.fTnI_PKA][end])
println("  fINa_PKA:   ", sol[sys.fINa_PKA][1], " → ", sol[sys.fINa_PKA][end])
println("  fINaK_PKA:  ", sol[sys.fINaK_PKA][1], " → ", sol[sys.fINaK_PKA][end])
println("  fRyR_PKA:   ", sol[sys.fRyR_PKA][1], " → ", sol[sys.fRyR_PKA][end])
println("  fIKur_PKA:  ", sol[sys.fIKur_PKA][1], " → ", sol[sys.fIKur_PKA][end])
