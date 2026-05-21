module GongBetaAdrenergicSignalingModelingToolkitExt

using GongBetaAdrenergicSignaling
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using PrecompileTools: @compile_workload

include("model.jl")

@compile_workload begin
    mtkcompile(GongBetaAdrenergic())
end

end
