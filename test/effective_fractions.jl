@testitem "Effective fractions (5 beats, iso_conc=0.1)" begin
    using GongBetaAdrenergicSignaling
    using ModelingToolkit
    using OrdinaryDiffEq
    using Test

    # Import the package function explicitly to avoid conflict with reference implementation
    import GongBetaAdrenergicSignaling: effective_fractions!

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
        3.1508784306208903e-01, 2.2764452443976033e-01, 1.1464122235575932e-02,
        3.1574761855345401e-01, 2.2830377222134693e-01, 1.2090971373466368e-02,
        6.5977549132490495e-04, 6.5924778156785065e-04, 6.2684913788924512e-04,
        8.0605594360421620e+00, 1.0168205888603646e+01, 7.2500064204073966e-01,
        1.5088670356824792e-02, 2.0301785223105639e-01, 9.4450105455938967e-03,
        3.4934216022719699e-05, 7.9258365871137903e-05, 4.6650531120865547e-06,
        1.5783661531335319e-01, 1.1198528787352134e-01, 6.5124400461221799e-01,
        4.5732893159637122e-01, 1.9391507301574226e-01, 1.9874885641697596e-01,
        1.7791124734365246e-01, 8.2450491583075913e-01, 5.7459191261599185e-01,
        2.4991300321677035e-01, 4.6468659270372151e-02, 7.2951454177341682e-02,
        5.0869426809509455e-01, 3.8174016338200079e-01, 1.2695410471288274e-01,
        3.4085311356084586e-02, 1.2856265289682841e-02, 9.1734505883606168e-03,
        4.2929159540311560e-05, 9.1818088807533484e-03, 1.2786328973410479e-02,
        1.3947070990702214e-03, 8.4982738942833240e-04, 6.7111102375698339e-01,
        6.7810850846464177e-01, 2.5472119057692644e-01, 1.4727233707191345e-01,
        7.7291029744852528e-03, 5.9774497788031483e-02, 2.7886391178291548e-02,
        1.3236022793859997e-04, 6.6134787186530375e-02, 6.6789462032286417e-02,
        6.5467484573846727e-04, 1.1024944246758660e-02, 6.2823949711590633e-06,
        1.8356368489705444e-02, 1.8995273934435657e-02, 6.3890544472593972e-04,
    ]

    # MATLAB reference data: phosphorylation fractions after 5 beats (8 values)
    matlab_phosph = [
        3.5035334571545548e-02,  # fICaLP
        5.1562755710366417e-03,  # fIKsP
        1.4958955220605257e-02,  # fPLBP
        1.4089877786878128e-02,  # fTnIP
        2.1446533563475978e-02,  # fINaP
        2.4002869024858656e-02,  # fINaKP
        3.4153168038786229e-02,  # fRyRP
        2.4984103654765344e-03,  # fIKurP
    ]

    # Build MTK model with isoproterenol
    @mtkcompile sys = GongBetaAdrenergic(iso_conc=iso_conc, radiusmultiplier=radiusmultiplier)

    # Create ODE problem and solve with Euler method
    prob_mtk = ODEProblem(sys, [], tspan)
    sol_mtk = solve(prob_mtk, Euler(); dt=1e-3)

    # Get reference initial state to create mapping
    u_ref = SignalingModel.get_default_initial_state()
    u0_mtk = prob_mtk.u0

    # Create mapping: ref_idx -> mtk_idx by matching initial conditions
    mapping = zeros(Int, length(u_ref))
    for i in 1:length(u_ref)
        for j in 1:length(u0_mtk)
            if isapprox(u_ref[i], u0_mtk[j], atol=1e-15)
                mapping[i] = j
                break
            end
        end
        @test mapping[i] != 0  # Ensure all states are mapped
    end

    # Extract final state in reference order
    julia_solution = sol_mtk.u[end][mapping]

    # Build reference parameters with same iso_conc/radiusmultiplier
    # Note: MTK parameter ordering differs from reference model c1-c167,
    # so we use reference parameters for effective_fractions! calculation
    c_ref = zeros(Float64, SignalingModel.NUM_PARAMS)
    SignalingModel.ConstantsSignalingMyokit2!(c_ref, iso_conc, radiusmultiplier)

    # Calculate effective fractions using package function
    effective_fractions = zeros(8)
    effective_fractions!(effective_fractions, julia_solution, c_ref)

    # Compare phosphorylation fractions to MATLAB reference
    @test all(isapprox.(effective_fractions, matlab_phosph, atol=1e-4))
end
