@testitem "Full integration (5 beats, iso_conc=0.1)" begin
    using GongBetaAdrenergicSignaling
    using ModelingToolkit
    using OrdinaryDiffEq
    using Test

    # Load reference model
    include("../references/julia/signaling.jl")
    using .SignalingModel

    # Test parameters - same as MATLAB version
    iso_conc = 0.1
    radiusmultiplier = 1.0
    beats = 5
    bcl = 1000
    tspan = (0.0, bcl * beats)

    # MATLAB reference data: final solution after 5 beats (57 state values)
    matlab_solution = [
        3.1508784306208903e-1,
        2.2764452443976033e-1,
        1.1464122235575932e-2,
        3.1574761855345401e-1,
        2.2830377222134693e-1,
        1.2090971373466368e-2,
        6.5977549132490495e-4,
        6.5924778156785065e-4,
        6.2684913788924512e-4,
        8.060559436042162e+0,
        1.0168205888603646e+1,
        7.2500064204073966e-1,
        1.5088670356824792e-2,
        2.0301785223105639e-1,
        9.4450105455938967e-3,
        3.4934216022719699e-5,
        7.9258365871137903e-5,
        4.6650531120865547e-6,
        1.5783661531335319e-1,
        1.1198528787352134e-1,
        6.5124400461221799e-1,
        4.5732893159637122e-1,
        1.9391507301574226e-1,
        1.9874885641697596e-1,
        1.7791124734365246e-1,
        8.2450491583075913e-1,
        5.7459191261599185e-1,
        2.4991300321677035e-1,
        4.6468659270372151e-2,
        7.2951454177341682e-2,
        5.0869426809509455e-1,
        3.8174016338200079e-1,
        1.2695410471288274e-1,
        3.4085311356084586e-2,
        1.2856265289682841e-2,
        9.1734505883606168e-3,
        4.292915954031156e-5,
        9.1818088807533484e-3,
        1.2786328973410479e-2,
        1.3947070990702214e-3,
        8.498273894283324e-4,
        6.7111102375698339e-1,
        6.7810850846464177e-1,
        2.5472119057692644e-1,
        1.4727233707191345e-1,
        7.7291029744852528e-3,
        5.9774497788031483e-2,
        2.7886391178291548e-2,
        1.3236022793859997e-4,
        6.6134787186530375e-2,
        6.6789462032286417e-2,
        6.5467484573846727e-4,
        1.102494424675866e-2,
        6.2823949711590633e-6,
        1.8356368489705444e-2,
        1.8995273934435657e-2,
        6.3890544472593972e-4,
    ]

    # Build MTK model with isoproterenol
    @mtkcompile sys = GongBetaAdrenergic(
        iso_conc = iso_conc, radiusmultiplier = radiusmultiplier
    )

    # Create ODE problem and solve with Euler method
    prob_mtk = ODEProblem(sys, [], tspan)
    sol_mtk = solve(prob_mtk, Euler(); dt = 1.0e-3)

    # Get reference initial state to create mapping
    u_ref = SignalingModel.get_default_initial_state()
    u0_mtk = prob_mtk.u0

    # Create mapping: ref_idx -> mtk_idx by matching initial conditions
    mapping = zeros(Int, length(u_ref))
    for i in 1:length(u_ref)
        for j in 1:length(u0_mtk)
            if isapprox(u_ref[i], u0_mtk[j], atol = 1.0e-15)
                mapping[i] = j
                break
            end
        end
        @test mapping[i] != 0  # Ensure all states are mapped
    end

    # Extract final state in reference order
    julia_solution = sol_mtk.u[end][mapping]

    # Compare final solution to MATLAB reference
    @test all(isapprox.(julia_solution, matlab_solution, atol = 1.0e-8))
end
