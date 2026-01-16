@testitem "RHS evaluation" begin
    using GongBetaAdrenergicSignaling
    using ModelingToolkit
    using OrdinaryDiffEq
    using Test

    # Load reference model
    include("../references/julia/signaling.jl")
    using .SignalingModel

    # Initialize reference model constants (167 parameters)
    c = zeros(Float64, SignalingModel.NUM_PARAMS)
    iso_conc = 0.0
    radiusmultiplier = 1.0
    SignalingModel.ConstantsSignalingMyokit2!(c, iso_conc, radiusmultiplier)

    # Get reference initial state (57 states)
    u_ref = SignalingModel.get_default_initial_state()

    # Allocate derivative array for reference
    du_ref = zeros(Float64, SignalingModel.NUM_STATES)

    # Evaluate reference RHS at t=0
    t_eval = 0.0
    SignalingModel.rhs_signaling!(du_ref, u_ref, c, t_eval)

    # Build MTK model
    @mtkcompile sys = GongBetaAdrenergic()

    # Create ODE problem with default initial conditions
    prob_mtk = ODEProblem(sys, [], (0.0, 1.0))

    # Allocate derivative array for MTK model
    du_mtk = zeros(length(prob_mtk.u0))

    # Evaluate MTK RHS at t=0
    prob_mtk.f(du_mtk, prob_mtk.u0, prob_mtk.p, 0.0)

    # MTK may reorder states during structural_simplify
    # Create mapping: ref_idx -> mtk_idx by matching initial conditions
    u0_ref = u_ref
    u0_mtk = prob_mtk.u0
    mapping = zeros(Int, length(u0_ref))

    for i in 1:length(u0_ref)
        for j in 1:length(u0_mtk)
            if isapprox(u0_ref[i], u0_mtk[j], atol=1e-15)
                mapping[i] = j
                break
            end
        end
        @test mapping[i] != 0  # Ensure all states are mapped
    end

    # Test that initial conditions match
    @test isapprox(u_ref, u0_mtk[mapping])
    # Test that derivatives match
    @test isapprox(du_ref, du_mtk[mapping])
end

@testitem "RHS evaluation with isoproterenol" begin
    using GongBetaAdrenergicSignaling
    using ModelingToolkit
    using OrdinaryDiffEq
    using Test

    # Load reference model
    include("../references/julia/signaling.jl")
    using .SignalingModel

    # Test with 1.0 μM isoproterenol
    iso_conc = 1.0
    radiusmultiplier = 1.0

    # Initialize reference model constants with isoproterenol
    c = zeros(Float64, SignalingModel.NUM_PARAMS)
    SignalingModel.ConstantsSignalingMyokit2!(c, iso_conc, radiusmultiplier)

    # Get reference initial state (should be same as baseline for initial conditions)
    u_ref = SignalingModel.get_default_initial_state()

    # Allocate derivative array for reference
    du_ref = zeros(Float64, SignalingModel.NUM_STATES)

    # Evaluate reference RHS at t=0 with isoproterenol
    t_eval = 0.0
    SignalingModel.rhs_signaling!(du_ref, u_ref, c, t_eval)

    # Build MTK model with isoproterenol
    @mtkcompile sys = GongBetaAdrenergic(
        iso_conc=iso_conc, radiusmultiplier=radiusmultiplier
    )

    # Create ODE problem with default initial conditions
    prob_mtk = ODEProblem(sys, [], (0.0, 1.0))

    # Allocate derivative array for MTK model
    du_mtk = zeros(length(prob_mtk.u0))

    # Evaluate MTK RHS at t=0
    prob_mtk.f(du_mtk, prob_mtk.u0, prob_mtk.p, 0.0)

    # MTK may reorder states during structural_simplify
    # Create mapping: ref_idx -> mtk_idx by matching initial conditions
    u0_ref = u_ref
    u0_mtk = prob_mtk.u0
    mapping = zeros(Int, length(u0_ref))

    for i in 1:length(u0_ref)
        for j in 1:length(u0_mtk)
            if isapprox(u0_ref[i], u0_mtk[j], atol=1e-15)
                mapping[i] = j
                break
            end
        end
        @test mapping[i] != 0  # Ensure all states are mapped
    end

    # Test that initial conditions match
    @test isapprox(u_ref, u0_mtk[mapping])

    # Test that derivatives match with isoproterenol
    @test isapprox(du_ref, du_mtk[mapping])

    # Verify that iso_conc parameter is correctly set (c1 should equal iso_conc)
    # Extract parameter c1 from the MTK problem using symbolic access
    c1_mtk = prob_mtk.ps[sys.c1]
    @test isapprox(c1_mtk, iso_conc)
end
