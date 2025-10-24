"""
    GongBetaAdrenergicSignaling

A Julia package providing a ModelingToolkit.jl implementation of the Gong et al. (2020)
beta-adrenergic signaling model for cardiac myocytes.

Reference:
Quantitative analysis of variability in an integrated model of human ventricular
electrophysiology and β-adrenergic signaling by Jingqi Q.X. Gong, Monica E. Susilo,
Anna Sher, Cynthia J. Musante, Eric A. Sobie
DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009

# Features
- **57 state variables**: G-protein signaling, cAMP dynamics, PKA activation, PDE phosphorylation,
  PP1 inhibition, and substrate phosphorylation
- **167 parameters**: All parameters include default values from the reference implementation
- **Default initial conditions**: Steady-state values for baseline (no isoproterenol) conditions
- **Observable outputs**: Effective phosphorylation fractions for 8 cardiac substrates (ICaL, IKs, PLB, TnI, INa, INaK, RyR, IKur)

# Exports
- `GongBetaAdrenergic`: The MTK model component (use with `@mtkbuild`)

# Basic Usage
```julia
using GongBetaAdrenergicSignaling
using ModelingToolkit, OrdinaryDiffEq

# Create model with defaults (iso = 0.0)
@mtkcompile sys = GongBetaAdrenergic()

# Create ODE problem with default initial conditions
prob = ODEProblem(sys, [], (0.0, 1000.0))

# Solve with appropriate stiff solver
sol = solve(prob, Rodas5P())
```

# Overriding Parameters
```julia
# !!! not yet implemented !!!
# Create model with isoproterenol stimulation
@mtkcompile sys = GongBetaAdrenergic(iso = 1.0)  # 1 μM isoproterenol

# Access observables after solving
using Plots
plot(sol, idxs=[sys.fICaLP, sys.fPLBP, sys.fRyRP],
     labels=["ICaL" "PLB" "RyR"],
     xlabel="Time (ms)", ylabel="Phosphorylation Fraction")
```
"""
module GongBetaAdrenergicSignaling

using Reexport
@reexport using ModelingToolkit

include("model.jl")

export GongBetaAdrenergic

end
