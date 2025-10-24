# GongBetaAdrenergicSignaling

[![Build Status](https://github.com/DerangedIons/GongBetaAdrenergicSignaling.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DerangedIons/GongBetaAdrenergicSignaling.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/DerangedIons/GongBetaAdrenergicSignaling.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DerangedIons/GongBetaAdrenergicSignaling.jl)

A Julia package providing a ModelingToolkit.jl implementation of the Gong et al. (2020) beta-adrenergic signaling model for cardiac myocytes.

**Reference:**
Quantitative analysis of variability in an integrated model of human ventricular electrophysiology and β-adrenergic signaling by Jingqi Q.X. Gong, Monica E. Susilo, Anna Sher, Cynthia J. Musante, Eric A. Sobie
DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009

## Usage

```julia
using GongBetaAdrenergicSignaling
using OrdinaryDiffEq

# Create model with defaults (iso = 0.0)
@mtkcompile sys = GongBetaAdrenergic()

# Create ODE problem with default initial conditions
prob = ODEProblem(sys, [], (0.0, 1000.0))

# Solve with appropriate stiff solver
sol = solve(prob, Rodas5P())
```

The model includes 57 state variables covering G-protein signaling, cAMP dynamics, PKA activation, PDE phosphorylation, PP1 inhibition, and substrate phosphorylation. It provides observable outputs for effective phosphorylation fractions of 8 cardiac substrates (ICaL, IKs, PLB, TnI, INa, INaK, RyR, IKur).
