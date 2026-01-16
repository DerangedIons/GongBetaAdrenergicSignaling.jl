# Migrated to MTK v11 programmatic API
# Original: @mtkmodel macro from MTK v9/v10
# Reference: Gong et al. (2020) beta-adrenergic signaling model

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
    GongBetaAdrenergic(; iso_conc=0.0, radiusmultiplier=1.0, name=:GongBetaAdrenergic)

Build the Gong et al. beta-adrenergic signaling model as a ModelingToolkit System.

# Arguments
- `iso_conc::Real=0.0`: Isoproterenol concentration (μM)
- `radiusmultiplier::Real=1.0`: Cell radius scaling factor
- `name::Symbol=:GongBetaAdrenergic`: System name

# Returns
- `System`: Unsimplified system - use `@mtkcompile` to compile

# Example
```julia
using GongBetaAdrenergicSignaling
using OrdinaryDiffEq

@mtkcompile sys = GongBetaAdrenergic(iso_conc=1.0)
prob = ODEProblem(sys, [], (0.0, 1000.0))
sol = solve(prob, Rodas5P())
```
"""
function GongBetaAdrenergic(; iso_conc = 0.0, radiusmultiplier = 1.0, name = :GongBetaAdrenergic)
    # Compute parameter values from structural parameters
    p = compute_parameters(iso_conc, radiusmultiplier)

    # Declare symbolic parameters with computed defaults
    @parameters begin
        c1 = p.c1
        c2 = p.c2
        c3 = p.c3
        c4 = p.c4
        c5 = p.c5
        c6 = p.c6
        c7 = p.c7
        c8 = p.c8
        c9 = p.c9
        c10 = p.c10
        c11 = p.c11
        c12 = p.c12
        c13 = p.c13
        c14 = p.c14
        c15 = p.c15
        c16 = p.c16
        c17 = p.c17
        c18 = p.c18
        c19 = p.c19
        c20 = p.c20
        c21 = p.c21
        c22 = p.c22
        c23 = p.c23
        c24 = p.c24
        c25 = p.c25
        c26 = p.c26
        c27 = p.c27
        c28 = p.c28
        c29 = p.c29
        c30 = p.c30
        c31 = p.c31
        c32 = p.c32
        c33 = p.c33
        c34 = p.c34
        c35 = p.c35
        c36 = p.c36
        c37 = p.c37
        c38 = p.c38
        c39 = p.c39
        c40 = p.c40
        c41 = p.c41
        c42 = p.c42
        c43 = p.c43
        c44 = p.c44
        c45 = p.c45
        c46 = p.c46
        c47 = p.c47
        c48 = p.c48
        c49 = p.c49
        c50 = p.c50
        c51 = p.c51
        c52 = p.c52
        c53 = p.c53
        c54 = p.c54
        c55 = p.c55
        c56 = p.c56
        c57 = p.c57
        c58 = p.c58
        c59 = p.c59
        c60 = p.c60
        c61 = p.c61
        c62 = p.c62
        c63 = p.c63
        c64 = p.c64
        c65 = p.c65
        c66 = p.c66
        c67 = p.c67
        c68 = p.c68
        c69 = p.c69
        c70 = p.c70
        c71 = p.c71
        c72 = p.c72
        c73 = p.c73
        c74 = p.c74
        c75 = p.c75
        c76 = p.c76
        c77 = p.c77
        c78 = p.c78
        c79 = p.c79
        c80 = p.c80
        c81 = p.c81
        c82 = p.c82
        c83 = p.c83
        c84 = p.c84
        c85 = p.c85
        c86 = p.c86
        c87 = p.c87
        c88 = p.c88
        c89 = p.c89
        c90 = p.c90
        c91 = p.c91
        c92 = p.c92
        c93 = p.c93
        c94 = p.c94
        c95 = p.c95
        c96 = p.c96
        c97 = p.c97
        c98 = p.c98
        c99 = p.c99
        c100 = p.c100
        c101 = p.c101
        c102 = p.c102
        c103 = p.c103
        c104 = p.c104
        c105 = p.c105
        c106 = p.c106
        c107 = p.c107
        c108 = p.c108
        c109 = p.c109
        c110 = p.c110
        c111 = p.c111
        c112 = p.c112
        c113 = p.c113
        c114 = p.c114
        c115 = p.c115
        c116 = p.c116
        c117 = p.c117
        c118 = p.c118
        c119 = p.c119
        c120 = p.c120
        c121 = p.c121
        c122 = p.c122
        c123 = p.c123
        c124 = p.c124
        c125 = p.c125
        c126 = p.c126
        c127 = p.c127
        c128 = p.c128
        c129 = p.c129
        c130 = p.c130
        c131 = p.c131
        c132 = p.c132
        c133 = p.c133
        c134 = p.c134
        c135 = p.c135
        c136 = p.c136
        c137 = p.c137
        c138 = p.c138
        c139 = p.c139
        c140 = p.c140
        c141 = p.c141
        c142 = p.c142
        c143 = p.c143
        c144 = p.c144
        c145 = p.c145
        c146 = p.c146
        c147 = p.c147
        c148 = p.c148
        c149 = p.c149
        c150 = p.c150
        c151 = p.c151
        c152 = p.c152
        c153 = p.c153
        c154 = p.c154
        c155 = p.c155
        c156 = p.c156
        c157 = p.c157
        c158 = p.c158
        c159 = p.c159
        c160 = p.c160
        c161 = p.c161
        c162 = p.c162
        c163 = p.c163
        c164 = p.c164
        c165 = p.c165
        c166 = p.c166
        c167 = p.c167
    end

    # Declare state variables with initial conditions (57 states)
    @variables begin
        beta_cav_Gs_aGTP(t) = 0.00685533455220118
        beta_eca_Gs_aGTP(t) = 0.0184630401160325
        beta_cyt_Gs_aGTP(t) = 0.000731426797266862
        beta_cav_Gs_bg(t) = 0.0074626834509494
        beta_eca_Gs_bg(t) = 0.0191020788696145
        beta_cyt_Gs_bg(t) = 0.00115141961261304
        beta_cav_Gs_aGDP(t) = 0.000607348898749425
        beta_eca_Gs_aGDP(t) = 0.000639038753581265
        beta_cyt_Gs_aGDP(t) = 0.00041999281534611
        cAMP_cav(t) = 0.344257659177271
        cAMP_eca(t) = 9.62262190945345
        cAMP_cyt(t) = 0.474028267773051
        beta_cav_Rb1_pka_tot(t) = 0.0148474493496437
        beta_eca_Rb1_pka_tot(t) = 0.2030159857254
        beta_cyt_Rb1_pka_tot(t) = 0.00944459882156118
        beta_cav_Rb1_grk_tot(t) = 1.76916170568022e-10
        beta_eca_Rb1_grk_tot(t) = 8.36801924219387e-10
        beta_cyt_Rb1_grk_tot(t) = 5.01719476559362e-11
        pka_cav_ARC(t) = 0.0898875954193269
        pka_cav_A2RC(t) = 0.00272422687231002
        pka_cav_A2R(t) = 0.225041998219388
        pka_cav_C(t) = 0.032238122010232
        pka_cav_PKIC(t) = 0.192803876209063
        pka_eca_ARC(t) = 0.205457169881723
        pka_eca_A2RC(t) = 0.174050238959555
        pka_eca_A2R(t) = 0.817148066343192
        pka_eca_C(t) = 0.567236181979005
        pka_eca_PKIC(t) = 0.249911884364108
        pka_cyt_ARC(t) = 0.0646981828366729
        pka_cyt_A2RC(t) = 0.0664977613514265
        pka_cyt_A2R(t) = 0.489057632872032
        pka_cyt_C(t) = 0.362107101574369
        pka_cyt_PKIC(t) = 0.126950531297531
        PDE3_P_cav(t) = 0.02339549924788
        PDE3_P_cyt(t) = 0.0128401592747216
        PDE4_P_cav(t) = 0.00629647926927854
        PDE4_P_eca(t) = 4.29166115483698e-5
        PDE4_P_cyt(t) = 0.00917030613498568
        inhib1_p(t) = 0.0123536190564101
        ICaLp(t) = 0.000664158274421826
        IKsp(t) = 0.000765842738691197
        iup_f_plb(t) = 0.666165471397222
        f_tni(t) = 0.673477756497978
        ina_f_ina(t) = 0.236980176272067
        f_inak(t) = 0.124710628782511
        RyRp(t) = 0.00404925913372347
        f_ikur(t) = 0.0589106047787742
        beta_cav_Rb2_pka_tot(t) = 0.0274407333314977
        beta_cav_Rb2_grk_tot(t) = 6.32124571143896e-10
        beta_cav_Gi_aGTP(t) = 0.00159025206300466
        beta_cav_Gi_bg(t) = 0.00209267447556971
        beta_cav_Gi_aGDP(t) = 0.000502422412564797
        beta_eca_Rb2_pka_tot(t) = 0.0110248493074472
        beta_eca_Rb2_grk_tot(t) = 8.04005829146876e-11
        beta_eca_Gi_aGTP(t) = 0.000364313646402573
        beta_eca_Gi_bg(t) = 0.000705654325530757
        beta_eca_Gi_aGDP(t) = 0.000341340679127998
    end

    # Declare intermediate algebraic variables (no defaults - MTK will compute)
    @variables begin
        beta_cav_Rb1_np_tot(t)
        beta_cav_Rb2_np_tot(t)
        beta_cav_Gi_abg(t)
        beta_cav_Gs_abg(t)
        beta_cav_Gs_f_d(t)
        beta_cav_Gs_f_b(t)
        beta_cav_Gs_f_c(t)
        beta_cav_Gs_f_rr(t)
        beta_cav_Gs_f_yr(t)
        beta_cav_Gs_f_yi(t)
        beta_cav_Gs_f_mag(t)
        beta_cav_Gs_f_arg(t)
        beta_cav_Gs_f_x(t)
        beta_cav_Gs_f_r(t)
        beta_cav_Gs_f_i(t)
        beta_cav_Gs_f(t)
        beta_cav_Rb1_f(t)
        beta_cav_LRb1(t)
        beta_cav_LRb1Gs(t)
        beta_cav_Rb2_f(t)
        beta_cav_LRb2(t)
        beta_cav_Rb2Gs(t)
        beta_cav_Rb1Gs(t)
        beta_cav_LRb2Gs(t)
        beta_eca_Gs_abg(t)
        beta_eca_Rb2_np_tot(t)
        beta_eca_Rb1_np_tot(t)
        beta_eca_Gs_f_d(t)
        beta_eca_Gs_f_c(t)
        beta_eca_Gs_f_b(t)
        beta_eca_Gs_f_rr(t)
        beta_eca_Gs_f_yi(t)
        beta_eca_Gs_f_yr(t)
        beta_eca_Gs_f_mag(t)
        beta_eca_Gs_f_arg(t)
        beta_eca_Gs_f_x(t)
        beta_eca_Gs_f_i(t)
        beta_eca_Gs_f_r(t)
        beta_eca_Gs_f(t)
        beta_eca_Rb1_f(t)
        beta_eca_Rb2_f(t)
        beta_eca_LRb1(t)
        beta_eca_LRb2(t)
        beta_eca_LRb2Gs(t)
        beta_eca_LRb1Gs(t)
        beta_eca_Rb2Gs(t)
        beta_eca_Rb1Gs(t)
        beta_eca_RGs_tot(t)
        beta_eca_LRGs_tot(t)
        beta_eca_Gi_abg(t)
        beta_eca_Rb2_pka_f_c(t)
        beta_eca_Rb2_pka_f_b(t)
        beta_eca_Rb2_pka_f(t)
        beta_cyt_Gs_abg(t)
        beta_cyt_Rb1_np_tot(t)
        beta_cyt_Rb1_np_f_b(t)
        beta_cyt_Rb1_np_f_c(t)
        Rb1_np_f(t)
        beta_cyt_Gs_f(t)
        LRb1_np(t)
        LRb1Gs_np(t)
        Rb1Gs_np(t)
        beta_cav_RGs_tot(t)
        beta_cav_LRGs_tot(t)
        beta_cav_Rb2_pka_f_c(t)
        beta_cav_Rb2_pka_f_b(t)
        beta_cav_Rb2_pka_f(t)
        beta_cav_Gi_f(t)
        beta_cav_Rb2Gi(t)
        beta_cav_LRb2Gi(t)
        beta_eca_Gi_f(t)
        beta_eca_Rb2Gi(t)
        beta_eca_LRb2Gi(t)
        pka_cav_RCf(t)
        pka_eca_RCf(t)
        pka_eca_dcAMP(t)
        pka_cyt_RCf(t)
        pka_cav_dcAMP(t)
        pka_cyt_dcAMP(t)
        dcAMP_PDE2_cyt(t)
        dcAMP_PDE2_eca(t)
        dcAMP_PDE2_cav(t)
        dcAMP_PDE3_cav(t)
        dcAMP_PDE4_cyt(t)
        dcAMP_PDE4_eca(t)
        dcAMP_PDE4_cav(t)
        dcAMP_PDE3_cyt(t)
        camp_cAMP_cyt_pde(t)
        camp_cAMP_cyt_j1(t)
        camp_cAMP_cyt_j2(t)
        camp_cAMP_eca_pde(t)
        camp_cAMP_eca_j2(t)
        camp_cAMP_eca_j1(t)
        camp_cAMP_cav_pde(t)
        camp_cAMP_cav_j2(t)
        camp_cAMP_cav_j1(t)
        ac_kAC47_cyt_gsa(t)
        kAC47_cyt(t)
        ac_kAC56_cav_gsa(t)
        gsi(t)
        kAC56_cav(t)
        ac_kAC47_eca_gsa(t)
        kAC47_eca(t)
        ac_kAC56_cyt_gsa(t)
        kAC56_cyt(t)
        dcAMP_AC47_cyt(t)
        dcAMP_AC56_cyt(t)
        dcAMP_AC56_cav(t)
        dcAMP_AC47_eca(t)
        pp1_PP1f_cyt_sum(t)
        PP1f_cyt(t)
        di(t)
        iks_sig_IKsp_dif(t)
        akap_sig_RyRp_dif(t)
        akap_sig_ICaLp_dif(t)
    end

    # Declare observable variables (no defaults)
    @variables begin
        fICaL_PKA(t)
        fIKs_PKA(t)
        fPLB_PKA(t)
        fTnI_PKA(t)
        fINa_PKA(t)
        fINaK_PKA(t)
        fRyR_PKA(t)
        fIKur_PKA(t)
        fMyBPC_PKA(t)
        Whole_cell_PP1(t)
    end

    # Build equations array
    eqs = [
        # CAVEOLAR COMPARTMENT
        # Total concentration of non-phosphorylated B1AR in the caveolar subspace
        beta_cav_Rb1_np_tot ~ c87 - beta_cav_Rb1_pka_tot - beta_cav_Rb1_grk_tot,
        # Total concentration of non-phosphorylated B2AR in the caveolar subspace
        beta_cav_Rb2_np_tot ~ c89 - beta_cav_Rb2_pka_tot - beta_cav_Rb2_grk_tot,
        # Concentration of Gi holoenzyme in the caveolar subspace
        beta_cav_Gi_abg ~ c63 * c59 * c5 - beta_cav_Gi_aGTP - beta_cav_Gi_aGDP,
        # Concentration of Gs holoenzyme in the caveolar subspace
        beta_cav_Gs_abg ~ c61 * c58 * c5 - beta_cav_Gs_aGTP - beta_cav_Gs_aGDP,
        beta_cav_Gs_f_d ~ beta_cav_Gs_abg * c93 / c92,
        beta_cav_Gs_f_b ~
            (c94 + c91) / c92 + beta_cav_Rb1_np_tot + beta_cav_Rb2_np_tot - beta_cav_Gs_abg,
        beta_cav_Gs_f_c ~
            (
            c91 * (beta_cav_Rb1_np_tot - beta_cav_Gs_abg) +
                c94 * (beta_cav_Rb2_np_tot - beta_cav_Gs_abg) +
                c93
        ) / c92,
        beta_cav_Gs_f_rr ~ (
            -beta_cav_Gs_f_d / 27.0 * beta_cav_Gs_f_b^3.0 -
                beta_cav_Gs_f_b * beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_c / 108.0 +
                beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_d / 6.0 +
                beta_cav_Gs_f_c^3.0 / 27.0 +
                beta_cav_Gs_f_d * beta_cav_Gs_f_d / 4.0
        ),
        beta_cav_Gs_f_yr ~ (
            sqrt(max(0.0, beta_cav_Gs_f_rr)) +
                beta_cav_Gs_f_d / 2.0 +
                beta_cav_Gs_f_b * beta_cav_Gs_f_c / 6.0 - beta_cav_Gs_f_b^3.0 / 27.0
        ),
        beta_cav_Gs_f_yi ~ sqrt(max(0.0, -beta_cav_Gs_f_rr)),
        beta_cav_Gs_f_mag ~
            (
            beta_cav_Gs_f_yr * beta_cav_Gs_f_yr + beta_cav_Gs_f_yi * beta_cav_Gs_f_yi
        )^(1.0 / 6.0),
        beta_cav_Gs_f_arg ~ atan(beta_cav_Gs_f_yi / beta_cav_Gs_f_yr) / 3.0,
        beta_cav_Gs_f_x ~
            (beta_cav_Gs_f_c / 3.0 - beta_cav_Gs_f_b * beta_cav_Gs_f_b / 9.0) /
            (beta_cav_Gs_f_mag * beta_cav_Gs_f_mag),
        beta_cav_Gs_f_r ~
            beta_cav_Gs_f_mag * cos(beta_cav_Gs_f_arg) * (1.0 - beta_cav_Gs_f_x) -
            beta_cav_Gs_f_b / 3.0,
        beta_cav_Gs_f_i ~
            beta_cav_Gs_f_mag * sin(beta_cav_Gs_f_arg) * (1.0 + beta_cav_Gs_f_x),
        # Concentration of free Gs in the caveolar subspace
        beta_cav_Gs_f ~
            sqrt(beta_cav_Gs_f_r * beta_cav_Gs_f_r + beta_cav_Gs_f_i * beta_cav_Gs_f_i),
        # Concentration of free non-phosphorylated beta1AR in caveolar subspace
        beta_cav_Rb1_f ~
            beta_cav_Rb1_np_tot / (1.0 + c1 / c65 + beta_cav_Gs_f * (c67 + c1) / (c66 * c67)),
        # Concentration of non-phosphorylated Ligand / Receptor1 complexes in caveolar subspace
        beta_cav_LRb1 ~ c1 * beta_cav_Rb1_f / c65,
        # Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in caveolar subspace
        beta_cav_LRb1Gs ~ c1 * beta_cav_Rb1_f * beta_cav_Gs_f / (c66 * c67),
        # Concentration of free non-phosphorylated beta2AR in caveolar subspace
        beta_cav_Rb2_f ~
            beta_cav_Rb2_np_tot / (1.0 + c1 / c72 + beta_cav_Gs_f * (c69 + c1) / (c71 * c69)),
        # Concentration of non-phosphorylated Ligand / Receptor2 complexes in caveolar subspace
        beta_cav_LRb2 ~ c1 * beta_cav_Rb2_f / c72,
        # Concentration of non-phosphorylated Receptor2 / G-protein complexes in caveolar subspace
        beta_cav_Rb2Gs ~ beta_cav_Rb2_f * beta_cav_Gs_f / c71,
        # Concentration of non-phosphorylated Receptor1 / G-protein complexes in caveolar subspace
        beta_cav_Rb1Gs ~ beta_cav_Rb1_f * beta_cav_Gs_f / c66,
        # Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in caveolar subspace
        beta_cav_LRb2Gs ~ c1 * beta_cav_Rb2_f * beta_cav_Gs_f / (c71 * c69),
        # Concentration of total PKA-phosphorylated beta1 receptors
        D(beta_cav_Rb1_pka_tot) ~
            0.001 * (c84 * pka_cav_C * beta_cav_Rb1_np_tot - c85 * beta_cav_Rb1_pka_tot),
        # Concentration of total GRK-phosphorylated beta1 receptors
        D(beta_cav_Rb1_grk_tot) ~
            0.001 *
            (c83 * c86 * (beta_cav_LRb1 + beta_cav_LRb1Gs) - c82 * beta_cav_Rb1_grk_tot),
        # Concentration of total PKA-phosphorylated beta2 receptors
        D(beta_cav_Rb2_pka_tot) ~
            0.001 * (c84 * pka_cav_C * beta_cav_Rb2_np_tot - c85 * beta_cav_Rb2_pka_tot),
        # Concentration of total GRK-phosphorylated beta2 receptors
        D(beta_cav_Rb2_grk_tot) ~
            0.001 *
            (c83 * c86 * (beta_cav_LRb2 + beta_cav_LRb2Gs) - c82 * beta_cav_Rb2_grk_tot),

        # EXTRACAVEOLAR COMPARTMENT
        # Concentration of Gs holoenzyme in the extracaveolar space
        beta_eca_Gs_abg ~ c60 * c58 * c6 - beta_eca_Gs_aGTP - beta_eca_Gs_aGDP,
        # Total concentration of non-phosphorylated B2AR in the extracaveolar space
        beta_eca_Rb2_np_tot ~ c97 - beta_eca_Rb2_pka_tot - beta_eca_Rb2_grk_tot,
        # Total concentration of non-phosphorylated B1AR in the extracaveolar space
        beta_eca_Rb1_np_tot ~ c98 - beta_eca_Rb1_pka_tot - beta_eca_Rb1_grk_tot,
        beta_eca_Gs_f_d ~ beta_eca_Gs_abg * c101 / c102,
        beta_eca_Gs_f_c ~
            (
            c103 * (beta_eca_Rb1_np_tot - beta_eca_Gs_abg) +
                c100 * (beta_eca_Rb2_np_tot - beta_eca_Gs_abg) +
                c101
        ) / c102,
        beta_eca_Gs_f_b ~
            (c100 + c103) / c102 + beta_eca_Rb1_np_tot + beta_eca_Rb2_np_tot - beta_eca_Gs_abg,
        beta_eca_Gs_f_rr ~ (
            -beta_eca_Gs_f_d / 27.0 * beta_eca_Gs_f_b^3.0 -
                beta_eca_Gs_f_b * beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_c / 108.0 +
                beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_d / 6.0 +
                beta_eca_Gs_f_c^3.0 / 27.0 +
                beta_eca_Gs_f_d * beta_eca_Gs_f_d / 4.0
        ),
        beta_eca_Gs_f_yi ~ sqrt(max(0.0, -beta_eca_Gs_f_rr)),
        beta_eca_Gs_f_yr ~ (
            sqrt(max(0.0, beta_eca_Gs_f_rr)) +
                beta_eca_Gs_f_d / 2.0 +
                beta_eca_Gs_f_b * beta_eca_Gs_f_c / 6.0 - beta_eca_Gs_f_b^3.0 / 27.0
        ),
        beta_eca_Gs_f_mag ~
            (
            beta_eca_Gs_f_yr * beta_eca_Gs_f_yr + beta_eca_Gs_f_yi * beta_eca_Gs_f_yi
        )^(1.0 / 6.0),
        beta_eca_Gs_f_arg ~ atan(beta_eca_Gs_f_yi / beta_eca_Gs_f_yr) / 3.0,
        beta_eca_Gs_f_x ~
            (beta_eca_Gs_f_c / 3.0 - beta_eca_Gs_f_b * beta_eca_Gs_f_b / 9.0) /
            (beta_eca_Gs_f_mag * beta_eca_Gs_f_mag),
        beta_eca_Gs_f_i ~
            beta_eca_Gs_f_mag * sin(beta_eca_Gs_f_arg) * (1.0 + beta_eca_Gs_f_x),
        beta_eca_Gs_f_r ~
            beta_eca_Gs_f_mag * cos(beta_eca_Gs_f_arg) * (1.0 - beta_eca_Gs_f_x) -
            beta_eca_Gs_f_b / 3.0,
        # Concentration of free Gs in the extracaveolar space
        beta_eca_Gs_f ~
            sqrt(beta_eca_Gs_f_r * beta_eca_Gs_f_r + beta_eca_Gs_f_i * beta_eca_Gs_f_i),
        # Concentration of free non-phosphorylated beta1AR in the extracaveolar space
        beta_eca_Rb1_f ~
            beta_eca_Rb1_np_tot / (1.0 + c1 / c65 + beta_eca_Gs_f * (c67 + c1) / (c66 * c67)),
        # Concentration of free non-phosphorylated beta2AR in the extracaveolar space
        beta_eca_Rb2_f ~
            beta_eca_Rb2_np_tot / (1.0 + c1 / c72 + beta_eca_Gs_f * (c69 + c1) / (c71 * c69)),
        # Concentration of non-phosphorylated Ligand / Receptor1 complexes in the extracaveolar space
        beta_eca_LRb1 ~ c1 * beta_eca_Rb1_f / c65,
        # Concentration of non-phosphorylated Ligand / Receptor2 complexes in the extracaveolar space
        beta_eca_LRb2 ~ c1 * beta_eca_Rb2_f / c72,
        # Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in the extracaveolar space
        beta_eca_LRb2Gs ~ c1 * beta_eca_Rb2_f * beta_eca_Gs_f / (c71 * c69),
        # Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in the extracaveolar space
        beta_eca_LRb1Gs ~ c1 * beta_eca_Rb1_f * beta_eca_Gs_f / (c66 * c67),
        # Concentration of non-phosphorylated Receptor2 / G-protein complexes in the extracaveolar space
        beta_eca_Rb2Gs ~ beta_eca_Rb2_f * beta_eca_Gs_f / c71,
        # Concentration of non-phosphorylated Receptor1 / G-protein complexes in the extracaveolar space
        beta_eca_Rb1Gs ~ beta_eca_Rb1_f * beta_eca_Gs_f / c66,
        beta_eca_RGs_tot ~ beta_eca_Rb1Gs + c96 * beta_eca_Rb2Gs,
        beta_eca_LRGs_tot ~ beta_eca_LRb1Gs + c96 * beta_eca_LRb2Gs,
        # Concentration of Gi holoenzyme in the extracaveolar space
        beta_eca_Gi_abg ~ c64 * c59 * c6 - beta_eca_Gi_aGTP - beta_eca_Gi_aGDP,
        beta_eca_Rb2_pka_f_c ~ -beta_eca_Rb2_pka_tot * c70 * c73,
        beta_eca_Rb2_pka_f_b ~ (
            beta_eca_Gi_abg * (c1 + c73) - beta_eca_Rb2_pka_tot * (c73 + c1) +
                c70 * c73 * (1.0 + c1 / c68)
        ),
        # Concentration of free PKA-phosphorylated beta2AR in the extracaveolar space
        beta_eca_Rb2_pka_f ~
            (
            -beta_eca_Rb2_pka_f_b + sqrt(
                beta_eca_Rb2_pka_f_b * beta_eca_Rb2_pka_f_b -
                    4.0 * c99 * beta_eca_Rb2_pka_f_c,
            )
        ) / (2.0 * c99),
        # Concentration of total PKA-phosphorylated beta1 receptors in the extracaveolar space
        D(beta_eca_Rb1_pka_tot) ~
            0.001 * (c84 * pka_eca_C * beta_eca_Rb1_np_tot - c85 * beta_eca_Rb1_pka_tot),
        # Concentration of total GRK-phosphorylated beta1 receptors in the extracaveolar space
        D(beta_eca_Rb1_grk_tot) ~
            0.001 *
            (c83 * c95 * (beta_eca_LRb1 + beta_eca_LRb1Gs) - c82 * beta_eca_Rb1_grk_tot),
        # Concentration of total PKA-phosphorylated beta2 receptors in the extracaveolar space
        D(beta_eca_Rb2_pka_tot) ~
            0.001 * (c84 * pka_eca_C * beta_eca_Rb2_np_tot - c85 * beta_eca_Rb2_pka_tot),
        # Concentration of total GRK-phosphorylated beta2 receptors in the extracaveolar space
        D(beta_eca_Rb2_grk_tot) ~
            0.001 *
            (c83 * c95 * (beta_eca_LRb2 + beta_eca_LRb2Gs) - c82 * beta_eca_Rb2_grk_tot),

        # CYTOPLASM COMPARTMENT
        # Concentration of Gs holoenzyme in the cytoplasm
        beta_cyt_Gs_abg ~ c62 * c58 * c7 - beta_cyt_Gs_aGTP - beta_cyt_Gs_aGDP,
        # Total concentration of non-phosphorylated beta-1 AR in the cytoplasm
        beta_cyt_Rb1_np_tot ~ c104 - beta_cyt_Rb1_pka_tot - beta_cyt_Rb1_grk_tot,
        beta_cyt_Rb1_np_f_b ~ (
            beta_cyt_Gs_abg * (c67 + c1) - beta_cyt_Rb1_np_tot * (c67 + c1) +
                c66 * c67 * (1.0 + c1 / c65)
        ),
        beta_cyt_Rb1_np_f_c ~ -beta_cyt_Rb1_np_tot * c67 * c66,
        # Concentration of free non-phosphorylated beta-1 AR in the cytoplasm
        Rb1_np_f ~
            (
            -beta_cyt_Rb1_np_f_b + sqrt(
                beta_cyt_Rb1_np_f_b * beta_cyt_Rb1_np_f_b - 4.0 * c106 * beta_cyt_Rb1_np_f_c
            )
        ) / (2.0 * c106),
        # Concentration of free (non-complexed) Gs in the cytoplasm
        beta_cyt_Gs_f ~ beta_cyt_Gs_abg / (1.0 + Rb1_np_f / c66 * (1.0 + c1 / c67)),
        # Concentration of non-phosphorylated ligand / beta-1 AR complexes in the cytoplasm
        LRb1_np ~ c1 * Rb1_np_f / c65,
        # Concentration of non-phosphorylated ligand / beta-1 AR / Gs complexes in the cytoplasm
        LRb1Gs_np ~ c1 * Rb1_np_f * beta_cyt_Gs_f / (c66 * c67),
        # Concentration of non-phosphorylated beta-1 AR / Gs complexes in the cytoplasm
        Rb1Gs_np ~ beta_cyt_Gs_f * Rb1_np_f / c66,
        # Concentration of total PKA-phosphorylated receptors in the cytoplasm
        D(beta_cyt_Rb1_pka_tot) ~
            0.001 * (c84 * pka_cyt_C * beta_cyt_Rb1_np_tot - c85 * beta_cyt_Rb1_pka_tot),
        # Concentration of total GRK-phosphorylated receptors in the cytoplasm
        D(beta_cyt_Rb1_grk_tot) ~
            0.001 * (c83 * c105 * (LRb1_np + LRb1Gs_np) - c82 * beta_cyt_Rb1_grk_tot),

        # G-PROTEIN ACTIVATION
        beta_cav_RGs_tot ~ beta_cav_Rb1Gs + c88 * beta_cav_Rb2Gs,
        beta_cav_LRGs_tot ~ beta_cav_LRb1Gs + c88 * beta_cav_LRb2Gs,
        beta_cav_Rb2_pka_f_c ~ -beta_cav_Rb2_pka_tot * c70 * c73,
        beta_cav_Rb2_pka_f_b ~ (
            beta_cav_Gi_abg * (c1 + c73) - beta_cav_Rb2_pka_tot * (c73 + c1) +
                c70 * c73 * (1.0 + c1 / c68)
        ),
        # Concentration of free PKA-phosphorylated beta2AR in the caveolar subspace
        beta_cav_Rb2_pka_f ~
            (
            -beta_cav_Rb2_pka_f_b + sqrt(
                beta_cav_Rb2_pka_f_b * beta_cav_Rb2_pka_f_b -
                    4.0 * c90 * beta_cav_Rb2_pka_f_c,
            )
        ) / (2.0 * c90),
        # Concentration of free (non-complexed) Gi in the caveolar subspace
        beta_cav_Gi_f ~
            beta_cav_Gi_abg / (1.0 + beta_cav_Rb2_pka_f / c70 * (1.0 + c1 / c73)),
        # Concentration of PKA phosphorylated b2AR/Gi complexes in caveolar subspace
        beta_cav_Rb2Gi ~ beta_cav_Rb2_pka_f * beta_cav_Gi_f / c70,
        # Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in caveolar subspace
        beta_cav_LRb2Gi ~ beta_cav_Rb2Gi * c1 / c73,
        # Concentration of free (non-complexed) Gi in the extracaveolar space
        beta_eca_Gi_f ~
            beta_eca_Gi_abg / (1.0 + beta_eca_Rb2_pka_f / c70 * (1.0 + c1 / c73)),
        # Concentration of PKA phosphorylated b2AR/Gi complexes in extracaveolar space
        beta_eca_Rb2Gi ~ beta_eca_Rb2_pka_f * beta_eca_Gi_f / c70,
        # Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in extracaveolar space
        beta_eca_LRb2Gi ~ c1 / c73 * beta_eca_Rb2Gi,
        # Concentration of active Gs alpha subunit in caveolar subspace (tri-phosphate)
        D(beta_cav_Gs_aGTP) ~
            0.001 * (c75 * beta_cav_RGs_tot + c74 * beta_cav_LRGs_tot - c78 * beta_cav_Gs_aGTP),
        # Concentration of active Gi alpha subunit in caveolar subspace
        D(beta_cav_Gi_aGTP) ~
            0.001 * (c77 * beta_cav_Rb2Gi + c76 * beta_cav_LRb2Gi - c79 * beta_cav_Gi_aGTP),
        # Concentration of active Gs alpha subunit in the extracaveolar space (tri-phosphate)
        D(beta_eca_Gs_aGTP) ~
            0.001 * (c75 * beta_eca_RGs_tot + c74 * beta_eca_LRGs_tot - c78 * beta_eca_Gs_aGTP),
        # Concentration of active Gi alpha subunit in the extracaveolar space
        D(beta_eca_Gi_aGTP) ~
            0.001 * (c77 * beta_eca_Rb2Gi + c76 * beta_eca_LRb2Gi - c79 * beta_eca_Gi_aGTP),
        # Concentration of active Gs alpha subunit in cytoplasm
        D(beta_cyt_Gs_aGTP) ~
            0.001 * (c75 * Rb1Gs_np + c74 * LRb1Gs_np - c78 * beta_cyt_Gs_aGTP),
        # Concentration of active Gs beta-gamma subunit in caveolar subspace
        D(beta_cav_Gs_bg) ~
            0.001 * (
            c75 * beta_cav_RGs_tot + c74 * beta_cav_LRGs_tot -
                c80 * beta_cav_Gs_bg * beta_cav_Gs_aGDP
        ),
        # Concentration of active Gi beta-gamma subunit in caveolar subspace
        D(beta_cav_Gi_bg) ~
            0.001 * (
            c77 * beta_cav_Rb2Gi + c76 * beta_cav_LRb2Gi -
                c81 * beta_cav_Gi_bg * beta_cav_Gi_aGDP
        ),
        # Concentration of active Gs beta-gamma subunit in the extracaveolar space
        D(beta_eca_Gs_bg) ~
            0.001 * (
            c75 * beta_eca_RGs_tot + c74 * beta_eca_LRGs_tot -
                c80 * beta_eca_Gs_bg * beta_eca_Gs_aGDP
        ),
        # Concentration of active Gi beta-gamma subunit in the extracaveolar space
        D(beta_eca_Gi_bg) ~
            0.001 * (
            c77 * beta_eca_Rb2Gi + c76 * beta_eca_LRb2Gi -
                c81 * beta_eca_Gi_bg * beta_eca_Gi_aGDP
        ),
        # Concentration of active Gs beta-gamma subunit in cytoplasm
        D(beta_cyt_Gs_bg) ~
            0.001 *
            (c75 * Rb1Gs_np + c74 * LRb1Gs_np - c80 * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP),
        # Concentration of inactive Gs alpha subunit in caveolar subspace (di-phosphate)
        D(beta_cav_Gs_aGDP) ~
            0.001 * (c78 * beta_cav_Gs_aGTP - c80 * beta_cav_Gs_bg * beta_cav_Gs_aGDP),
        # Concentration of inactive Gi alpha subunit in caveolar subspace
        D(beta_cav_Gi_aGDP) ~
            0.001 * (c79 * beta_cav_Gi_aGTP - c81 * beta_cav_Gi_bg * beta_cav_Gi_aGDP),
        # Concentration of inactive Gs alpha subunit in the extracaveolar space (di-phosphate)
        D(beta_eca_Gs_aGDP) ~
            0.001 * (c78 * beta_eca_Gs_aGTP - c80 * beta_eca_Gs_bg * beta_eca_Gs_aGDP),
        # Concentration of inactive Gi alpha subunit in the extracaveolar space
        D(beta_eca_Gi_aGDP) ~
            0.001 * (c79 * beta_eca_Gi_aGTP - c81 * beta_eca_Gi_bg * beta_eca_Gi_aGDP),
        # Concentration of inactive Gs alpha subunit in cytoplasm
        D(beta_cyt_Gs_aGDP) ~
            0.001 * (c78 * beta_cyt_Gs_aGTP - c80 * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP),

        # PKA DYNAMICS
        # Concentration of free PKA RC subunits in the caveolar compartment
        pka_cav_RCf ~ c12 - pka_cav_ARC - pka_cav_A2RC - pka_cav_A2R,
        # Caveolar concentration of PKA RC dimer with 1 cAMP molecule bound
        D(pka_cav_ARC) ~
            0.001 * (
            c19 * pka_cav_RCf * cAMP_cav - c22 * pka_cav_ARC -
                c20 * pka_cav_ARC * cAMP_cav + c23 * pka_cav_A2RC
        ),
        # Caveolar concentration of PKA RC dimer with 2 cAMP molecules bound
        D(pka_cav_A2RC) ~
            0.001 * (
            c20 * pka_cav_ARC * cAMP_cav - (c23 + c21) * pka_cav_A2RC +
                c24 * pka_cav_A2R * pka_cav_C
        ),
        # Caveolar concentration of PKA R subunit with 2 cAMP molecules bound
        D(pka_cav_A2R) ~ 0.001 * (c21 * pka_cav_A2RC - c24 * pka_cav_A2R * pka_cav_C),
        # Caveolar concentration of free PKA catalytic subunit
        D(pka_cav_C) ~
            0.001 * (
            c21 * pka_cav_A2RC - c24 * pka_cav_A2R * pka_cav_C + c18 * pka_cav_PKIC -
                c17 * (c14 - pka_cav_PKIC) * pka_cav_C
        ),
        # Caveolar concentration of free PKI inactivated PKA C subunit
        D(pka_cav_PKIC) ~
            0.001 * (c17 * (c14 - pka_cav_PKIC) * pka_cav_C - c18 * pka_cav_PKIC),
        # Concentration of free PKA RC subunits in the Extracaveolar compartment
        pka_eca_RCf ~ c11 - pka_eca_ARC - pka_eca_A2RC - pka_eca_A2R,
        # Extracaveolar rate of change in free cAMP through binding by PKA
        pka_eca_dcAMP ~ (
            -c19 * pka_eca_RCf * cAMP_eca + c25 * pka_eca_ARC -
                c20 * pka_eca_ARC * cAMP_eca + c26 * pka_eca_A2RC
        ),
        # Extracaveolar concentration of PKA RC dimer with 1 cAMP molecule bound
        D(pka_eca_ARC) ~
            0.001 * (
            c19 * pka_eca_RCf * cAMP_eca - c25 * pka_eca_ARC -
                c20 * pka_eca_ARC * cAMP_eca + c26 * pka_eca_A2RC
        ),
        # Extracaveolar concentration of PKA RC dimer with 2 cAMP molecules bound
        D(pka_eca_A2RC) ~
            0.001 * (
            c20 * pka_eca_ARC * cAMP_eca - (c26 + c21) * pka_eca_A2RC +
                c27 * pka_eca_A2R * pka_eca_C
        ),
        # Extracaveolar concentration of PKA R subunit with 2 cAMP molecules bound
        D(pka_eca_A2R) ~ 0.001 * (c21 * pka_eca_A2RC - c27 * pka_eca_A2R * pka_eca_C),
        # Extracaveolar concentration of free PKA catalytic subunit
        D(pka_eca_C) ~
            0.001 * (
            c21 * pka_eca_A2RC - c27 * pka_eca_A2R * pka_eca_C + c18 * pka_eca_PKIC -
                c17 * (c15 - pka_eca_PKIC) * pka_eca_C
        ),
        # Extracaveolar concentration of free PKI inactivated PKA C subunit
        D(pka_eca_PKIC) ~
            0.001 * (c17 * (c15 - pka_eca_PKIC) * pka_eca_C - c18 * pka_eca_PKIC),
        # Concentration of free PKA RC subunits in the Cytosolic compartment
        pka_cyt_RCf ~ c13 - pka_cyt_ARC - pka_cyt_A2RC - pka_cyt_A2R,
        # Cytosolic concentration of PKA RC dimer with 1 cAMP molecule bound
        D(pka_cyt_ARC) ~
            0.001 * (
            c19 * pka_cyt_RCf * cAMP_cyt - c28 * pka_cyt_ARC -
                c20 * pka_cyt_ARC * cAMP_cyt + c29 * pka_cyt_A2RC
        ),
        # Cytosolic concentration of PKA RC dimer with 2 cAMP molecules bound
        D(pka_cyt_A2RC) ~
            0.001 * (
            c20 * pka_cyt_ARC * cAMP_cyt - (c29 + c21) * pka_cyt_A2RC +
                c30 * pka_cyt_A2R * pka_cyt_C
        ),
        # Cytosolic concentration of PKA R subunit with 2 cAMP molecules bound
        D(pka_cyt_A2R) ~ 0.001 * (c21 * pka_cyt_A2RC - c30 * pka_cyt_A2R * pka_cyt_C),
        # Cytosolic concentration of free PKA catalytic subunit
        D(pka_cyt_C) ~
            0.001 * (
            c21 * pka_cyt_A2RC - c30 * pka_cyt_A2R * pka_cyt_C + c18 * pka_cyt_PKIC -
                c17 * (c16 - pka_cyt_PKIC) * pka_cyt_C
        ),
        # Cytosolic concentration of free PKI inactivated PKA C subunit
        D(pka_cyt_PKIC) ~
            0.001 * (c17 * (c16 - pka_cyt_PKIC) * pka_cyt_C - c18 * pka_cyt_PKIC),

        # cAMP DYNAMICS
        # Caveolar rate of change in free cAMP through binding by PKA
        pka_cav_dcAMP ~ (
            -c19 * pka_cav_RCf * cAMP_cav + c22 * pka_cav_ARC -
                c20 * pka_cav_ARC * cAMP_cav + c23 * pka_cav_A2RC
        ),
        # Cytosolic rate of change in free cAMP through binding by PKA
        pka_cyt_dcAMP ~ (
            -c19 * pka_cyt_RCf * cAMP_cyt + c28 * pka_cyt_ARC -
                c20 * pka_cyt_ARC * cAMP_cyt + c29 * pka_cyt_A2RC
        ),

        # PDE DEGRADATION
        # Rate of cAMP degradation by PDE2 in cytosolic subspace
        dcAMP_PDE2_cyt ~ c52 * c41 / (1.0 + c44 / cAMP_cyt),
        # Rate of cAMP degradation by PDE2 in extracaveolar subspace
        dcAMP_PDE2_eca ~ c51 * c41 / (1.0 + c44 / cAMP_eca),
        # Rate of cAMP degradation by PDE2 in caveolar subspace
        dcAMP_PDE2_cav ~ c50 * c41 / (1.0 + c44 / cAMP_cav),
        # Rate of cAMP degradation by PDE3 in caveolar subspace
        dcAMP_PDE3_cav ~ (c53 + (c47 - 1.0) * PDE3_P_cav) * c42 / (1.0 + c45 / cAMP_cav),
        # Rate of cAMP degradation by PDE4 in cytosolic subspace
        dcAMP_PDE4_cyt ~ (c57 + (c47 - 1.0) * PDE4_P_cyt) * c43 / (1.0 + c46 / cAMP_cyt),
        # Rate of cAMP degradation by PDE4 in extracaveolar subspace
        dcAMP_PDE4_eca ~ (c56 + (c47 - 1.0) * PDE4_P_eca) * c43 / (1.0 + c46 / cAMP_eca),
        # Rate of cAMP degradation by PDE4 in caveolar subspace
        dcAMP_PDE4_cav ~ (c55 + (c47 - 1.0) * PDE4_P_cav) * c43 / (1.0 + c46 / cAMP_cav),
        # Rate of cAMP degradation by PDE3 in cytosolic subspace
        dcAMP_PDE3_cyt ~ (c54 + (c47 - 1.0) * PDE3_P_cyt) * c42 / (1.0 + c45 / cAMP_cyt),
        camp_cAMP_cyt_pde ~ dcAMP_PDE2_cyt + dcAMP_PDE3_cyt + dcAMP_PDE4_cyt,
        camp_cAMP_cyt_j1 ~ c9 * (cAMP_cav - cAMP_cyt) / c4,
        camp_cAMP_cyt_j2 ~ c10 * (cAMP_eca - cAMP_cyt) / c4,
        camp_cAMP_eca_pde ~ dcAMP_PDE2_eca + dcAMP_PDE4_eca,
        camp_cAMP_eca_j2 ~ c10 * (cAMP_eca - cAMP_cyt) / c3,
        camp_cAMP_eca_j1 ~ c8 * (cAMP_cav - cAMP_eca) / c3,
        camp_cAMP_cav_pde ~ dcAMP_PDE2_cav + dcAMP_PDE3_cav + dcAMP_PDE4_cav,
        camp_cAMP_cav_j2 ~ c9 * (cAMP_cav - cAMP_cyt) / c2,
        camp_cAMP_cav_j1 ~ c8 * (cAMP_cav - cAMP_eca) / c2,

        # ADENYLYL CYCLASE
        ac_kAC47_cyt_gsa ~ beta_cyt_Gs_aGTP^c107,
        kAC47_cyt ~ c116 * (c114 + ac_kAC47_cyt_gsa / (c110 + ac_kAC47_cyt_gsa)),
        ac_kAC56_cav_gsa ~ beta_cav_Gs_aGTP^c108,
        gsi ~ beta_cav_Gs_aGTP^c109,
        kAC56_cav ~ (
            c117 *
                (c115 + ac_kAC56_cav_gsa / (c111 + ac_kAC56_cav_gsa)) *
                (
                1.0 -
                    (1.0 - c118 * gsi / (c113 + gsi)) * beta_cav_Gi_bg /
                    (c112 + beta_cav_Gi_bg)
            )
        ),
        ac_kAC47_eca_gsa ~ beta_eca_Gs_aGTP^c107,
        kAC47_eca ~ c116 * (c114 + ac_kAC47_eca_gsa / (c110 + ac_kAC47_eca_gsa)),
        ac_kAC56_cyt_gsa ~ beta_cyt_Gs_aGTP^c108,
        kAC56_cyt ~ c117 * (c115 + ac_kAC56_cyt_gsa / (c111 + ac_kAC56_cyt_gsa)),
        # Rate of cAMP production by AC type 4/7 in cytoplasm
        dcAMP_AC47_cyt ~ kAC47_cyt * c119 * c121,
        # Rate of cAMP production by AC type 5/6 in cytoplasm
        dcAMP_AC56_cyt ~ kAC56_cyt * c123 * c121,
        # Rate of cAMP production by AC type 5/6 in caveolar subspace
        dcAMP_AC56_cav ~ kAC56_cav * c120 * c121,
        # Rate of cAMP production by AC type 4/7 in extracaveolar subspace
        dcAMP_AC47_eca ~ kAC47_eca * c122 * c121,
        # Caveolar concentration of cAMP
        D(cAMP_cav) ~
            0.001 * (
            pka_cav_dcAMP + dcAMP_AC56_cav - camp_cAMP_cav_pde - camp_cAMP_cav_j1 -
                camp_cAMP_cav_j2
        ),
        # Extracaveolar concentration of cAMP
        D(cAMP_eca) ~
            0.001 * (
            pka_eca_dcAMP + dcAMP_AC47_eca - camp_cAMP_eca_pde + camp_cAMP_eca_j1 -
                camp_cAMP_eca_j2
        ),
        # Cytosolic concentration of cAMP
        D(cAMP_cyt) ~
            0.001 * (
            pka_cyt_dcAMP + dcAMP_AC47_cyt + dcAMP_AC56_cyt - camp_cAMP_cyt_pde +
                camp_cAMP_cyt_j1 +
                camp_cAMP_cyt_j2
        ),

        # PDE PHOSPHORYLATION
        # Concentration of phosphorylated PDE3 in the caveolar subspace
        D(PDE3_P_cav) ~ 0.001 * (c48 * pka_cav_C * (c53 - PDE3_P_cav) - c49 * PDE3_P_cav),
        # Concentration of phosphorylated PDE3 in the cytosolic subspace
        D(PDE3_P_cyt) ~ 0.001 * (c48 * pka_cyt_C * (c54 - PDE3_P_cyt) - c49 * PDE3_P_cyt),
        # Concentration of phosphorylated PDE4 in the caveolar subspace
        D(PDE4_P_cav) ~ 0.001 * (c48 * pka_cav_C * (c55 - PDE4_P_cav) - c49 * PDE4_P_cav),
        # Concentration of phosphorylated PDE4 in the extracaveolar subspace
        D(PDE4_P_eca) ~ 0.001 * (c48 * pka_eca_C * (c56 - PDE4_P_eca) - c49 * PDE4_P_eca),
        # Concentration of phosphorylated PDE4 in the cytosolic subspace
        D(PDE4_P_cyt) ~ 0.001 * (c48 * pka_cyt_C * (c57 - PDE4_P_cyt) - c49 * PDE4_P_cyt),

        # PP1 INHIBITION
        pp1_PP1f_cyt_sum ~ c38 - c37 + inhib1_p,
        # Concentration of uninhibited PP1 in the cytosolic compartment
        PP1f_cyt ~ 0.5 * (sqrt(pp1_PP1f_cyt_sum^2.0 + 4.0 * c38 * c37) - pp1_PP1f_cyt_sum),
        di ~ c39 - inhib1_p,
        # Concentration of phosphorylated PP1 inhibitor 1 (cytoplasmic)
        D(inhib1_p) ~
            0.001 *
            (c31 * pka_cyt_C * di / (c33 + di) - c32 * c40 * inhib1_p / (c34 + inhib1_p)),

        # CHANNEL PHOSPHORYLATION - Substrates without AKAP
        # Fraction of phosphorylated PLB
        D(iup_f_plb) ~
            0.001 * (
            c132 * pka_cyt_C * (1.0 - iup_f_plb) / (c134 + 1.0 - iup_f_plb) -
                c133 * PP1f_cyt * iup_f_plb / (c135 + iup_f_plb)
        ),
        # Fraction of phosphorylated Troponin
        D(f_tni) ~
            0.001 * (
            c147 * pka_cyt_C * (1.0 - f_tni) / (c149 + 1.0 - f_tni) -
                c148 * c40 * f_tni / (c150 + f_tni)
        ),
        # Fraction of phosphorylated INa channels
        D(ina_f_ina) ~
            0.001 * (
            c130 * pka_cav_C * (1.0 - ina_f_ina) / (c128 + 1.0 - ina_f_ina) -
                c131 * c36 * ina_f_ina / (c129 + ina_f_ina)
        ),
        # Fraction of phosphorylated INaK
        D(f_inak) ~
            0.001 * (
            c124 * pka_cav_C * (1.0 - f_inak) / (c126 + 1.0 - f_inak) -
                c125 * c36 * f_inak / (c127 + f_inak)
        ),
        # Fraction of phosphorylated IKur channels
        D(f_ikur) ~
            0.001 * (
            c136 * pka_eca_C * (1.0 - f_ikur) / (c138 + 1.0 - f_ikur) -
                c137 * c35 * f_ikur / (c139 + f_ikur)
        ),

        # CHANNEL PHOSPHORYLATION - Substrates with AKAP
        iks_sig_IKsp_dif ~ c146 - IKsp,
        # Concentration of phosphorylated IKs channels
        D(IKsp) ~
            0.001 * (
            c140 * pka_eca_C * iks_sig_IKsp_dif / (c142 + iks_sig_IKsp_dif) -
                c141 * c35 * IKsp / (c143 + IKsp)
        ),
        akap_sig_RyRp_dif ~ c162 - RyRp,
        D(RyRp) ~
            0.001 * (
            c152 * pka_cav_C * akap_sig_RyRp_dif / (c154 + akap_sig_RyRp_dif) -
                c153 * c36 * RyRp / (c155 + RyRp)
        ),
        akap_sig_ICaLp_dif ~ c164 - ICaLp,
        # Concentration of phosphorylated L-type Calcium channels
        D(ICaLp) ~
            0.001 * (
            c157 * pka_cav_C * akap_sig_ICaLp_dif / (c159 + akap_sig_ICaLp_dif) -
                c158 * c36 * ICaLp / (c160 + ICaLp)
        ),

        # EFFECTIVE PHOSPHORYLATION FRACTIONS (Observables)
        # ICaL effective fraction (L-type calcium channel)
        fICaL_PKA ~ let
            fp = (ICaLp + c163) / c156
            fp_clamped = ifelse(fp < 0.0001, 0.0001, ifelse(fp > 0.9999, 0.9999, fp))
            raw = (fp_clamped - c166) / (0.9273 - c166)
            ifelse(raw < 0.0, 0.0, ifelse(raw > 1.0, 1.0, raw))
        end,

        # IKs effective fraction (slow delayed rectifier K+ channel)
        fIKs_PKA ~ let
            fp = (IKsp + c145) / c144
            fp_clamped = ifelse(fp < 0.0001, 0.0001, ifelse(fp > 0.9999, 0.9999, fp))
            raw = (fp_clamped - c167) / (0.785 - c167)
            ifelse(raw < 0.0, 0.0, ifelse(raw > 1.0, 1.0, raw))
        end,

        # PLB effective fraction (phospholamban, affects SERCA)
        fPLB_PKA ~ let
            raw = (iup_f_plb - 0.6662) / (0.9945 - 0.6662)
            ifelse(raw < 0.0, 0.0, ifelse(raw > 1.0, 1.0, raw))
        end,

        # TnI effective fraction (troponin I, affects calcium sensitivity)
        fTnI_PKA ~ let
            raw = (f_tni - 0.6735188) / (0.9991797 - 0.6735188)
            ifelse(raw < 0.0, 0.0, ifelse(raw > 1.0, 1.0, raw))
        end,

        # INa effective fraction (sodium channel)
        fINa_PKA ~ let
            raw = (ina_f_ina - 0.2394795) / (0.9501431 - 0.2394795)
            ifelse(raw < 0.0, 0.0, ifelse(raw > 1.0, 1.0, raw))
        end,

        # INaK effective fraction (Na+/K+ pump)
        fINaK_PKA ~ let
            raw = (f_inak - 0.1263453) / (0.9980137 - 0.1263453)
            ifelse(raw < 0.0, 0.0, ifelse(raw > 1.0, 1.0, raw))
        end,

        # RyR effective fraction (ryanodine receptor, calcium release)
        fRyR_PKA ~ let
            fp = (RyRp + c161) / c151
            fp_clamped = ifelse(fp < 0.0001, 0.0001, ifelse(fp > 0.9999, 0.9999, fp))
            raw = (fp_clamped - c165) / (0.9586 - c165)
            ifelse(raw < 0.0, 0.0, ifelse(raw > 1.0, 1.0, raw))
        end,

        # IKur effective fraction (ultra-rapid delayed rectifier K+ channel)
        fIKur_PKA ~ let
            raw = (f_ikur - 5.893798e-2) / (0.393747 - 5.893798e-2)
            ifelse(raw < 0.0, 0.0, ifelse(raw > 1.0, 1.0, raw))
        end,

        # MyBPC effective fraction (myosin binding protein C, uses same as TnI)
        fMyBPC_PKA ~ fTnI_PKA,

        # Whole cell PP1 concentration (protein phosphatase 1)
        Whole_cell_PP1 ~ c36 / c5 + c35 / c6 + PP1f_cyt / c7,
    ]

    # Return the System
    return System(eqs, t; name)
end
