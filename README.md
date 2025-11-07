# GongBetaAdrenergicSignaling

[![Build Status](https://github.com/DerangedIons/GongBetaAdrenergicSignaling.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DerangedIons/GongBetaAdrenergicSignaling.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package providing a ModelingToolkit.jl implementation of the Gong et al. (2020) beta-adrenergic signaling model for cardiac myocytes.

## Features

- **57 state variables**: G-protein signaling, cAMP dynamics, PKA activation, PDE phosphorylation, PP1 inhibition, and substrate phosphorylation
- **167 parameters**: Automatically computed from structural parameters
- **Isoproterenol stimulation**: Simply set `iso_conc` and all 167 parameters automatically update
- **Phosphorylation observables**: Built-in observables for effective fractions of 8 cardiac ion channels/proteins
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

## Phosphorylation Fraction Observables

The model automatically computes effective phosphorylation fractions for 8 cardiac substrates as **observables**. These are available directly from the solution:

```julia
using GongBetaAdrenergicSignaling
using OrdinaryDiffEq

# Solve the model
@mtkcompile sys = GongBetaAdrenergic(iso_conc=1.0)
prob = ODEProblem(sys, [], (0.0, 5000.0))
sol = solve(prob, Rodas5P())

# Access phosphorylation fractions (available at all time points)
fICaL_PKA = sol[sys.fICaL_PKA]   # L-type calcium channel
fIKs_PKA = sol[sys.fIKs_PKA]     # Slow delayed rectifier K+ channel
fPLB_PKA = sol[sys.fPLB_PKA]     # Phospholamban (SERCA regulation)
fTnI_PKA = sol[sys.fTnI_PKA]     # Troponin I (myofilament Ca2+ sensitivity)
fINa_PKA = sol[sys.fINa_PKA]     # Sodium channel
fINaK_PKA = sol[sys.fINaK_PKA]   # Na+/K+ pump
fRyR_PKA = sol[sys.fRyR_PKA]     # Ryanodine receptor (Ca2+ release)
fIKur_PKA = sol[sys.fIKur_PKA]   # Ultra-rapid delayed rectifier K+ channel
fMyBPC_PKA = sol[sys.fMyBPC_PKA] # Myosin binding protein C

# Extract final values
final_fICaL = sol[sys.fICaL_PKA][end]
```

These observables are computed automatically during the solve and represent the effective fraction of each substrate that is phosphorylated by PKA.

## Reference

Gong, J.Q.X., Susilo, M.E., Sher, A., Musante, C.J., & Sobie, E.A. (2020).
Quantitative analysis of variability in an integrated model of human ventricular
electrophysiology and β-adrenergic signaling. *Journal of Molecular and Cellular Cardiology*, 143, 96-106.
DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009
