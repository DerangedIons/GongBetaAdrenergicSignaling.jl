"""
    GongBetaAdrenergicSignaling

Plain-Julia implementation of the Gong et al. (2020) beta-adrenergic signaling model for
cardiac myocytes (57 states, 167 parameters, 8 phosphorylation observables).

The model is exposed as an in-place ODE right-hand side, [`rhs_signaling!`](@ref), for
direct use with OrdinaryDiffEq.jl:

```julia
using GongBetaAdrenergicSignaling, OrdinaryDiffEq

p = compute_parameters(1.0)                  # 1 μM isoproterenol
u0 = default_initial_state()
prob = ODEProblem(rhs_signaling!, u0, (0.0, 1000.0), p)
sol = solve(prob, Rodas5P())
```

A symbolic ModelingToolkit `System` is also available through [`GongBetaAdrenergic`](@ref),
which lives in a package extension — load it with `using ModelingToolkit`.

# Reference
Gong, J.Q.X., Susilo, M.E., Sher, A., Musante, C.J., & Sobie, E.A. (2020).
Quantitative analysis of variability in an integrated model of human ventricular
electrophysiology and β-adrenergic signaling. Journal of Molecular and Cellular Cardiology.
DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009
"""
module GongBetaAdrenergicSignaling

using PrecompileTools: @compile_workload

"Number of differential state variables in the signaling model."
const NUM_STATES = 57

"Number of parameters (c1-c167) in the signaling model."
const NUM_PARAMS = 167

include("parameters.jl")
include("rhs.jl")
include("effective_fractions.jl")
include("initial_conditions.jl")

"""
    GongBetaAdrenergic(; iso_conc=0.0, radiusmultiplier=1.0, name=:GongBetaAdrenergic)

Build the Gong et al. beta-adrenergic signaling model as a ModelingToolkit `System`.

!!! note
    This constructor is provided by a package extension. Run `using ModelingToolkit`
    before calling it, otherwise a `MethodError` is thrown. For the plain ODE function
    that needs no extra dependencies, use [`rhs_signaling!`](@ref) instead.

# Arguments
- `iso_conc=0.0`: Isoproterenol concentration (μM)
- `radiusmultiplier=1.0`: Cell radius scaling factor
- `name=:GongBetaAdrenergic`: System name

# Returns
- `System`: unsimplified system — pass it through `mtkcompile` before solving.
"""
function GongBetaAdrenergic end

export rhs_signaling!,
    effective_fractions!,
    compute_parameters,
    compute_parameters!,
    default_initial_state,
    NUM_STATES,
    NUM_PARAMS,
    GongBetaAdrenergic

@compile_workload begin
    p = compute_parameters(1.0)
    u0 = default_initial_state()
    du = similar(u0)
    rhs_signaling!(du, u0, p, 0.0)
    effective_fractions!(zeros(8), u0, p)
end

end
