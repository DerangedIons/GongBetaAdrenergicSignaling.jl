"""
Basic simulation example for the Gong Beta Adrenergic Signaling model.

This example demonstrates:
1. Baseline simulation (no isoproterenol)
2. Isoproterenol stimulation (1.0 μM)
3. Visualization of key phosphorylation targets
"""

using GongBetaAdrenergicSignaling
using ModelingToolkit
using OrdinaryDiffEq
using CairoMakie

N_beats = 200
BCL = 1000

# Baseline simulation
@mtkbuild sys = GongBetaAdrenergic()
prob = ODEProblem(sys, [], (0.0, N_beats * BCL), [])
sol = solve(prob, Tsit5())

## Plotting
fig = Figure(size=(600, 400));

ax1 = Axis(fig[1, 1],
           xlabel="Time (ms)",
           ylabel="Phosphorylation Fraction",
           title="Baseline (no ISO)")
#lines!(ax1, sol.t, sol[sys.ICaLp], label="ICaL", linewidth=2)
#lines!(ax1, sol.t, sol[sys.IKsp], label="IKs", linewidth=2)
#lines!(ax1, sol.t, sol[sys.RyRp], label="RyR", linewidth=2)
lines!(ax1, sol.t, sol[sys.ICaLp], label="ICaL", linewidth=2)
#lines!(ax1, sol.t, sol[sys.iup_f_plb], label="PLB", linewidth=2)
axislegend(ax1, position=:rb)

display(fig)
