"""
    GongBetaAdrenergicSignaling

ModelingToolkit.jl implementation of the Gong et al. (2020) beta-adrenergic signaling
model for cardiac myocytes. The model contains 57 states, 167 parameters, and computes
effective phosphorylation fractions for 8 cardiac ion channels and regulatory proteins.

# Exports
- `GongBetaAdrenergic`: MTK model constructor with built-in observable outputs for 8 phosphorylation fractions
- `compute_parameters`: Compute parameters from `iso_conc` and `radiusmultiplier`

See the README for usage examples and documentation.

# Reference
Gong, J.Q.X., Susilo, M.E., Sher, A., Musante, C.J., & Sobie, E.A. (2020).
Quantitative analysis of variability in an integrated model of human ventricular
electrophysiology and β-adrenergic signaling. Journal of Molecular and Cellular Cardiology.
DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009
"""
module GongBetaAdrenergicSignaling

using Reexport
@reexport using ModelingToolkit

# Load parameter computation function first
include("parameters.jl")

# Load the @mtkmodel definition (exports GongBetaAdrenergic directly)
include("model.jl")

export GongBetaAdrenergic, compute_parameters

# Precompilation workload
using PrecompileTools: @compile_workload

@compile_workload begin
    # Precompile parameter computation (baseline and with iso)
    compute_parameters()
    compute_parameters(1.0, 1.0)

    # Precompile model construction
    sys_baseline = GongBetaAdrenergic()
    compiled_baseline = mtkcompile(sys_baseline)

    sys_stim = GongBetaAdrenergic(iso_conc=1.0)
    compiled_stim = mtkcompile(sys_stim)
end

end
