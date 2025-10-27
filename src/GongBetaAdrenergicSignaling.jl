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
- **167 parameters**: Automatically computed from structural parameters `iso_conc` and `radiusmultiplier`
- **Default initial conditions**: Steady-state values for baseline (no isoproterenol) conditions
- **Observable outputs**: Effective phosphorylation fractions for 8 cardiac substrates (ICaL, IKs, PLB, TnI, INa, INaK, RyR, IKur)

# Exports
- `GongBetaAdrenergic`: MTK model with structural parameters `iso_conc` (μM) and `radiusmultiplier`
- `compute_parameters`: Standalone utility to compute all 167 parameters from `iso_conc` and `radiusmultiplier`

# Basic Usage
```julia
using GongBetaAdrenergicSignaling
using OrdinaryDiffEq

# Create model with defaults (iso_conc=0.0, radiusmultiplier=1.0)
@mtkcompile sys = GongBetaAdrenergic()

# Create ODE problem with default initial conditions
prob = ODEProblem(sys, [], (0.0, 1000.0))

# Solve with appropriate stiff solver
sol = solve(prob, Rodas5P())
```

# With Isoproterenol Stimulation
```julia
# Use structural parameters to set iso_conc
@mtkcompile sys = GongBetaAdrenergic(iso_conc=1.0)

# Or customize both structural parameters
@mtkcompile sys = GongBetaAdrenergic(iso_conc=1.0, radiusmultiplier=1.2)

# Create ODE problem and solve
prob = ODEProblem(sys, [], (0.0, 1000.0))
sol = solve(prob, Rodas5P())
```

# Advanced: Manual Parameter Computation
```julia
# Compute parameters explicitly (for inspection or external use)
params = compute_parameters(1.0, 1.0)

# Inspect specific parameter (e.g., c1 = iso_conc)
params.c1  # Returns 1.0

# Note: Parameters are now computed internally via structural parameters,
# but you can still use compute_parameters() as a standalone utility
```
"""
module GongBetaAdrenergicSignaling

using Reexport
@reexport using ModelingToolkit

# Load parameter computation function first
include("parameters.jl")

# Load the @mtkmodel definition (exports GongBetaAdrenergic directly)
include("model.jl")

export GongBetaAdrenergic, compute_parameters

end
