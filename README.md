# GongBetaAdrenergicSignaling

[![Build Status](https://github.com/DerangedIons/GongBetaAdrenergicSignaling.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DerangedIons/GongBetaAdrenergicSignaling.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/DerangedIons/GongBetaAdrenergicSignaling.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DerangedIons/GongBetaAdrenergicSignaling.jl)

A Julia package providing a ModelingToolkit.jl implementation of the Gong et al. (2020) beta-adrenergic signaling model for cardiac myocytes.

## Features

- **57 state variables**: G-protein signaling, cAMP dynamics, PKA activation, PDE phosphorylation, PP1 inhibition, and substrate phosphorylation
- **167 parameters**: Automatically computed from structural parameters
- **Isoproterenol stimulation**: Simply set `iso_conc` and all 167 parameters automatically update
- **Phosphorylation outputs**: Compute effective fractions for 8 cardiac ion channels/proteins
- **ModelingToolkit integration**: Symbolic model construction with automatic code generation

## Basic Usage

```julia
using GongBetaAdrenergicSignaling
using OrdinaryDiffEq

# Create model with baseline (no isoproterenol)
@mtkcompile sys = GongBetaAdrenergic()

# Create ODE problem with default initial conditions
prob = ODEProblem(sys, [], (0.0, 1000.0))

# Solve with stiff solver
sol = solve(prob, Rodas5P())
```

## Isoproterenol Stimulation

The key feature of this model is the ability to simulate beta-adrenergic stimulation by simply changing the `iso_conc` parameter. **All 167 model parameters automatically recalculate** based on the isoproterenol concentration:

```julia
# Simulate with 1.0 μM isoproterenol
@mtkcompile sys = GongBetaAdrenergic(iso_conc=1.0)
prob = ODEProblem(sys, [], (0.0, 1000.0))
sol = solve(prob, Rodas5P())

# Try different concentrations
for iso in [0.0, 0.01, 0.1, 1.0]
    @mtkcompile sys = GongBetaAdrenergic(iso_conc=iso)
    prob = ODEProblem(sys, [], (0.0, 5000.0))
    sol = solve(prob, Rodas5P())
    # Analyze dose-response...
end
```

You can also adjust the cell geometry via `radiusmultiplier`:

```julia
@mtkcompile sys = GongBetaAdrenergic(iso_conc=1.0, radiusmultiplier=1.2)
```

## Computing Phosphorylation Fractions

The model outputs can be post-processed to compute effective phosphorylation fractions for cardiac substrates:

```julia
using GongBetaAdrenergicSignaling

# After solving the model...
fractions = zeros(8)
effective_fractions!(fractions, states, parameters)

# fractions now contains:
# [1] ICaL  - L-type calcium channel
# [2] IKs   - Slow delayed rectifier K+ channel
# [3] PLB   - Phospholamban (SERCA regulation)
# [4] TnI   - Troponin I (myofilament Ca2+ sensitivity)
# [5] INa   - Sodium channel
# [6] INaK  - Na+/K+ pump
# [7] RyR   - Ryanodine receptor (Ca2+ release)
# [8] IKur  - Ultra-rapid delayed rectifier K+ channel
```

## Reference

Gong, J.Q.X., Susilo, M.E., Sher, A., Musante, C.J., & Sobie, E.A. (2020).
Quantitative analysis of variability in an integrated model of human ventricular
electrophysiology and β-adrenergic signaling. *Journal of Molecular and Cellular Cardiology*, 143, 96-106.
DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009
