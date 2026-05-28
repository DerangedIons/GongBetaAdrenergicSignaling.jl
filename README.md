# GongBetaAdrenergicSignaling

[![Build Status](https://github.com/DerangedIons/GongBetaAdrenergicSignaling.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DerangedIons/GongBetaAdrenergicSignaling.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia implementation of the Gong et al. (2020) beta-adrenergic signaling model for cardiac myocytes.

## Features

- **57 state variables**: G-protein signaling, cAMP dynamics, PKA activation, PDE phosphorylation, PP1 inhibition, and substrate phosphorylation
- **167 parameters**: Automatically computed from structural parameters
- **Isoproterenol stimulation**: Simply set `iso_conc` and all 167 parameters automatically update
- **Phosphorylation fractions**: Effective fractions for 8 cardiac ion channels/proteins
- **Plain ODE function**: A standard in-place `f(du, u, p, t)` right-hand side — no heavy dependencies, solved directly by OrdinaryDiffEq.jl
- **Optional ModelingToolkit model**: A symbolic `System` is available through a package extension

## Basic Usage

The model is exposed as [`rhs_signaling!`](src/rhs.jl), an in-place ODE right-hand side
that plugs straight into OrdinaryDiffEq.jl:

```julia
using GongBetaAdrenergicSignaling
using OrdinaryDiffEq

p = compute_parameters()                 # 167-element parameter vector (baseline)
u0 = default_initial_state()             # 57-element initial state

prob = ODEProblem(rhs_signaling!, u0, (0.0, 1000.0), p)
sol = solve(prob, Rodas5P())             # stiff solver recommended
```

## Isoproterenol Stimulation

Beta-adrenergic stimulation is set through the `iso_conc` argument of `compute_parameters`.
**All 167 model parameters automatically recalculate** based on the isoproterenol concentration:

```julia
# Simulate with 1.0 μM isoproterenol
p = compute_parameters(1.0)
prob = ODEProblem(rhs_signaling!, default_initial_state(), (0.0, 1000.0), p)
sol = solve(prob, Rodas5P())

# Dose-response sweep
for iso in (0.0, 0.01, 0.1, 1.0)
    p = compute_parameters(iso)
    prob = ODEProblem(rhs_signaling!, default_initial_state(), (0.0, 5000.0), p)
    sol = solve(prob, Rodas5P())
    # analyze...
end
```

Cell geometry can be adjusted via the second argument, `radiusmultiplier`:

```julia
p = compute_parameters(1.0, 1.2)
```

For repeated parameter builds, `compute_parameters!(c, iso_conc, radiusmultiplier)` fills a
preallocated 167-element vector in place.

## Phosphorylation Fractions

`effective_fractions!` computes the effective phosphorylation fractions of 8 cardiac
substrates from a state vector:

```julia
sol = solve(ODEProblem(rhs_signaling!, default_initial_state(), (0.0, 5000.0),
                       compute_parameters(1.0)), Rodas5P())

# Order: [fICaL, fIKs, fPLB, fTnI, fINa, fINaK, fRyR, fIKur]
out = zeros(8)
effective_fractions!(out, sol.u[end], compute_parameters(1.0))
```

| Index | Substrate | Description |
|-------|-----------|-------------|
| 1 | ICaL  | L-type calcium channel |
| 2 | IKs   | Slow delayed rectifier K⁺ channel |
| 3 | PLB   | Phospholamban (SERCA regulation) |
| 4 | TnI   | Troponin I (myofilament Ca²⁺ sensitivity) |
| 5 | INa   | Sodium channel |
| 6 | INaK  | Na⁺/K⁺ pump |
| 7 | RyR   | Ryanodine receptor (Ca²⁺ release) |
| 8 | IKur  | Ultra-rapid delayed rectifier K⁺ channel |

## ModelingToolkit Model (optional)

A symbolic ModelingToolkit `System` is available through `GongBetaAdrenergic`. It lives in a
package extension, so ModelingToolkit is a *weak* dependency — it is not loaded (and does not
precompile) unless you ask for it:

```julia
using GongBetaAdrenergicSignaling
using ModelingToolkit          # activates the extension
using OrdinaryDiffEq

sys = mtkcompile(GongBetaAdrenergic(iso_conc = 1.0))
prob = ODEProblem(sys, [], (0.0, 1000.0))
sol = solve(prob, Rodas5P())
```

Calling `GongBetaAdrenergic` without `using ModelingToolkit` throws a `MethodError`.

## Reference

Gong, J.Q.X., Susilo, M.E., Sher, A., Musante, C.J., & Sobie, E.A. (2020).
Quantitative analysis of variability in an integrated model of human ventricular
electrophysiology and β-adrenergic signaling. *Journal of Molecular and Cellular Cardiology*, 143, 96-106.
DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009
