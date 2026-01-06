@testitem "Effective fractions (5 beats, iso_conc=0.1)" begin
    using GongBetaAdrenergicSignaling
    using ModelingToolkit
    using OrdinaryDiffEq
    using Test

    # Test parameters - same as MATLAB version
    iso_conc = 0.1
    radiusmultiplier = 1.0
    beats = 5
    bcl = 1000
    tspan = (0.0, bcl * beats)

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
    @mtkcompile sys = GongBetaAdrenergic(
        iso_conc=iso_conc, radiusmultiplier=radiusmultiplier
    )

    # Create ODE problem and solve with Euler method
    prob_mtk = ODEProblem(sys, [], tspan)
    sol_mtk = solve(prob_mtk, Euler(); dt=1e-3)

    # Extract effective fractions from model observables at final time
    effective_fractions = [
        sol_mtk[sys.fICaL_PKA][end],
        sol_mtk[sys.fIKs_PKA][end],
        sol_mtk[sys.fPLB_PKA][end],
        sol_mtk[sys.fTnI_PKA][end],
        sol_mtk[sys.fINa_PKA][end],
        sol_mtk[sys.fINaK_PKA][end],
        sol_mtk[sys.fRyR_PKA][end],
        sol_mtk[sys.fIKur_PKA][end],
    ]

    # Compare phosphorylation fractions to MATLAB reference
    @test all(isapprox.(effective_fractions, matlab_phosph, atol=1e-4))
end
