@testitem "Effective fractions (5 beats, iso_conc=0.1)" begin
    using GongBetaAdrenergicSignaling
    using OrdinaryDiffEqLowOrderRK
    using Test

    iso_conc = 0.1
    radiusmultiplier = 1.0
    beats = 5
    bcl = 1000
    tspan = (0.0, bcl * beats)

    # MATLAB reference data: phosphorylation fractions after 5 beats (8 values)
    matlab_phosph = [
        3.5035334571545548e-2,  # fICaLP
        5.1562755710366417e-3,  # fIKsP
        1.4958955220605257e-2,  # fPLBP
        1.4089877786878128e-2,  # fTnIP
        2.1446533563475978e-2,  # fINaP
        2.4002869024858656e-2,  # fINaKP
        3.4153168038786229e-2,  # fRyRP
        2.4984103654765344e-3,  # fIKurP
    ]

    # Solve the plain ODE function over 5 beats with the Euler method
    p = compute_parameters(iso_conc, radiusmultiplier)
    u0 = default_initial_state()
    prob = ODEProblem(rhs_signaling!, u0, tspan, p)
    sol = solve(prob, Euler(); dt = 1.0e-3)

    # Compute effective phosphorylation fractions at the final state
    out = zeros(8)
    effective_fractions!(out, sol.u[end], p)

    @test all(isapprox.(out, matlab_phosph, atol = 1.0e-4))
end
