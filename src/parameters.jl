"""
    compute_parameters(iso_conc=0.0, radiusmultiplier=1.0)

Compute all 167 model parameters (c1-c167) based on isoproterenol concentration and radius multiplier.

This function implements the same parameter calculation logic as `ConstantsSignalingMyokit2!`
from the reference model, translating drug concentrations and geometric properties into
model-specific rate constants, equilibrium constants, and concentrations.

# Arguments
- `iso_conc::Real=0.0`: Isoproterenol concentration (μM)
- `radiusmultiplier::Real=1.0`: Cell radius scaling factor

# Returns
- `NamedTuple`: Parameters c1 through c167 as named tuple for passing to model constructor

# Example
```julia
# Default baseline (no stimulation)
params = compute_parameters()

# 1 μM isoproterenol stimulation
params = compute_parameters(1.0, 1.0)

# Use with model constructor
@mtkcompile sys = GongBetaAdrenergic(; params...)
```
"""
function compute_parameters(iso_conc::Real=0.0, radiusmultiplier::Real=1.0)
    c = zeros(Float64, 167)

    # iso
    c[1] = iso_conc
    # IBMX
    ibmx = 0.0  # Concentration of IBMX (in the range 0 to 100) not used in OharaBA

    # Cell geometry
    length = 0.01  # Cell length
    cell_pi = 3.14159265358979312e00  # pi
    radius = 0.0011 * radiusmultiplier  # Cell radius
    volume = 1000.0 * cell_pi * radius * radius * length  # Cell volume

    c[2] = 0.02 * volume  # Volume of the caveolar subspace
    c[3] = 0.04 * volume  # Volume of the extracaveolar subspace
    c[4] = volume * 0.678  # Volume of the Cytoplasm / Myoplasm

    c[5] = volume / c[2]  # Ratio of whole volume to caveolar subspace volume
    c[6] = volume / c[3]  # Ratio of whole volume to extracaveolar subspace volume
    c[7] = volume / c[4]  # Ratio of whole volume to cytoplasm volume

    # cAMP
    c[8] = 5e-15 * 1000000.0  # Rate of cAMP diffusion between caveolar and cytosolic compartments
    c[9] = 7.5e-14 * 1000000.0  # Rate of cAMP diffusion between caveolar and extracaveolar compartments
    c[10] = 9e-15 * 1000000.0  # Rate of cAMP diffusion between extracaveolar and cytosolic compartments

    # PKA (See pg 29 - 32 of Heijman supplementary document)
    PKA_tot = 0.5  # Total cellular concentration of PKA holoenzyme (umol/L)
    f_cav = 0.0388  # Fraction of PKA located in caveolar compartment
    f_eca = 0.1  # Fraction of PKA located in extracaveolar compartment
    f_cyt = 1.0 - f_cav - f_eca  # Fraction of PKA located in cytosolic compartment

    c[11] = f_eca * PKA_tot * c[6]  # Concentration of PKA in the extracaveolar compartment
    c[12] = f_cav * PKA_tot * c[5]  # Concentration of PKA in the caveolar compartment
    c[13] = f_cyt * PKA_tot * c[7]  # Concentration of PKA in the cytosolic compartment

    PKI_tot = 0.2 * PKA_tot  # Total cellular concentration of PKA inhibitor
    f_pki_cav = f_cav  # Fraction of PKI located in caveolar compartment
    f_pki_eca = f_eca  # Fraction of PKI located in extracaveolar compartment
    f_pki_cyt = 1.0 - f_pki_cav - f_pki_eca  # Fraction of PKI located in cytosolic compartment

    c[14] = f_pki_cav * PKI_tot * c[5]  # Concentration of protein kinase inhibitor in caveolar compartment
    c[15] = f_pki_eca * PKI_tot * c[6]  # Concentration of protein kinase inhibitor in extracaveolar compartment
    c[16] = f_pki_cyt * PKI_tot * c[7]  # Concentration of protein kinase inhibitor in cytosolic compartment

    c[17] = 50.0  # Forward rate for inhibition of C subunit by PKI
    K_pki = 0.01 / 50.0  # Equilibrium value for inhibition of C subunit by PKI (umol/L)
    c[18] = c[17] * K_pki  # Backward rate for inhibition of C subunit by PKI

    # pka_cav
    c[19] = 100.0  # Caveolar forward rate for binding of the first cAMP to PKA
    c[20] = 100.0  # Caveolar forward rate for binding of the second cAMP to PKA
    c[21] = 100.0  # Caveolar forward rate for dissociation of C subunit

    pka_cav_K1 = 2.4984  # Caveolar equilibrium value for the binding of the first cAMP to PKA
    pka_cav_K2 = 11.359  # Caveolar equilibrium value for the binding of the second cAMP to PKA
    pka_cav_K3 = 0.3755  # Caveolar equilibrium value for dissociation of C subunit

    c[22] = c[19] * pka_cav_K1  # Caveolar backward rate for binding of the first cAMP to PKA
    c[23] = c[20] * pka_cav_K2  # Caveolar backward rate for binding of the second cAMP to PKA
    c[24] = c[21] * pka_cav_K3  # Caveolar backward rate for dissociation of C subunit

    pka_eca_K1 = pka_cav_K1  # Extracaveolar equilibrium value for the binding of the first cAMP to PKA
    pka_eca_K2 = pka_cav_K2  # Extracaveolar equilibrium value for the binding of the second cAMP to PKA
    pka_eca_K3 = pka_cav_K3  # Extracaveolar equilibrium value for dissociation of C subunit

    c[25] = c[19] * pka_eca_K1  # Extracaveolar backward rate for binding of the first cAMP to PKA
    c[26] = c[20] * pka_eca_K2  # Extracaveolar backward rate for binding of the second cAMP to PKA
    c[27] = c[21] * pka_eca_K3  # Extracaveolar backward rate for dissociation of C subunit

    pka_cyt_K1 = 0.1088  # Cytosolic equilibrium value for the binding of the first cAMP to PKA
    pka_cyt_K2 = 0.4612  # Cytosolic equilibrium value for binding of the second cAMP to PKA
    pka_cyt_K3 = 0.3755  # Cytosolic equilibrium value for dissociation of C subunit

    c[28] = c[19] * pka_cyt_K1  # Cytosolic backward rate for binding of the first cAMP to PKA
    c[29] = c[20] * pka_cyt_K2  # Cytosolic backward rate for binding of the second cAMP to PKA
    c[30] = c[21] * pka_cyt_K3  # Cytosolic backward rate for dissociation of C subunit

    # PP1 (see pg 32-33 of Heijman Supplementary document)
    f = 0.3  # Fractional increase in PP1 after Inh1 knockout
    c[31] = 0.010145  # Rate of phosphorylation of inhibitor 1 by PKA
    c[32] = 0.0035731  # Rate of dephosphorylation if inhibitor 1
    c[33] = 0.001469  # Affinity of inhibitor 1 for PKA catalytic subunit
    c[34] = 1.95259999999999991e-05  # Affinity of inhibitor 1 for PP2A

    c[35] = 0.1  # PP1 concentration in the extracaveolar compartment
    c[36] = 0.25  # PP1 concentration in the caveolar compartment
    c[37] = 0.2  # PP1 concentration in the cytosolic compartment

    c[38] = 0.001  # Affinity foor PP1 Inhibitor 1 binding

    c[39] = f / (1.0 - f) * c[38] + f * c[37]  # Concentration of phosphatase inhibitor 1 in the cytosolic compartment
    c[40] = 1.0  # PP2A concentration?

    # PDE (see pg 26 - 29 of Heijman supplementary document)
    PDE2_tot = 0.029268  # Total cellular concentration of PDE2
    f_pde2_cav = 0.16957  # Fraction of PDE2 located in caveolar compartment
    f_pde2_eca = 2.12570000000000006e-04  # Fraction of PDE2 located in extracaveolar compartment
    f_pde2_cyt = 1.0 - f_pde2_cav - f_pde2_eca  # Fraction of PDE2 located in cytosolic compartment
    f_pde2_part = f_pde2_cav + f_pde2_eca  # Fraction of PDE2 located in the "particulate fraction" (cav + eca)

    f_pde_part = 0.2  # Fraction of total PDE located in the "particulate fraction" (cav + eca)
    f_pde4_part = 0.125  # Fraction of PDE4 located in the "particulate fraction" (cav + eca)

    f_pde4_cav = 0.12481  # Fraction of PDE4 in caveolar compartment
    f_pde4_eca = f_pde4_part - f_pde4_cav  # Fraction of PDE4 located in the extracaveolar compartment
    f_pde4_cyt = 1.0 - f_pde4_part  # Fraction of PDE4 in cytosolic compartment

    c[41] = 20.0  # Rate of cAMP hydrolysis by PDE2
    c[42] = 2.5  # Rate of cAMP hydrolysis by PDE3
    c[43] = 4.0  # Rate of cAMP hydrolysis by PDE4

    c[44] = 50.0  # Affinity of PDE2 for cAMP
    c[45] = 0.8  # Affinity of PDE3 for cAMP
    c[46] = 1.4  # Affinity of PDE4 for cAMP

    KmIbmxPde2 = 21.58  # Affinity of IBMX for PDE2
    h_ibmx_pde2 = 1.167  # Hill coefficient for inhibition of PDE2 by IBMX
    h_ibmx_pde3 = 0.7629  # Hill coefficient for inhibition of PDE3 by IBMX
    h_ibmx_pde4 = 0.9024  # Hill coefficient for inhibition of PDE4 by IBMX

    KmIbmxPde3 = 2.642  # Affinity of IBMX for PDE3
    KmIbmxPde4 = 11.89  # Affinity of IBMX for PDE4

    KPDEp = 0.52218
    c[47] = 3.0  # Increase in PDE3 / PDE4 activity after phosphorylation
    ff_pde3_cyt = 0.35  # Fraction of PDE in cytosol that is of type 3
    c[48] = 0.0196  # Rate of phosphorylation by PKA of PDE3 and PDE4
    r_pde34_frac = 3.71  # Ratio of PDE3 to PDE4 in particulate fraction (78:21)
    r_pde3_cyt = ff_pde3_cyt / (1.0 - ff_pde3_cyt)  # Relative contribution of PDE3 to cytosolic PDE3
    c[49] = KPDEp * c[48]  # Rate of dephosphorylation of PDE3 and PDE4

    pde_PDE3_tot_alpha =
        r_pde3_cyt * (
            f_pde4_part * (1.0 + r_pde34_frac - r_pde34_frac * f_pde2_part - f_pde_part) +
            f_pde2_part * (f_pde_part - 1.0)
        ) + r_pde34_frac * f_pde4_part * (f_pde_part - f_pde2_part)
    pde_PDE3_tot_beta =
        f_pde4_part * (1.0 + r_pde34_frac + f_pde_part * (r_pde3_cyt - r_pde34_frac)) -
        f_pde_part * (1.0 + r_pde3_cyt)

    PDE3_tot = pde_PDE3_tot_alpha / pde_PDE3_tot_beta * PDE2_tot  # Total cellular concentration of PDE3
    PDE4_tot =
        ((f_pde_part - f_pde2_part) * PDE2_tot + f_pde_part * PDE3_tot) /
        ((1.0 + r_pde34_frac) * f_pde4_part - f_pde_part)  # Total cellular concentration of PDE4

    ibmx_h2 = ibmx^h_ibmx_pde2
    ibmx_h3 = ibmx^h_ibmx_pde3
    ibmx_h4 = ibmx^h_ibmx_pde4

    ibmx2 = (1.0 - ibmx_h2 / (KmIbmxPde2 + ibmx_h2)) * PDE2_tot
    ibmx3 = (1.0 - ibmx_h3 / (KmIbmxPde3 + ibmx_h3)) * PDE3_tot
    ibmx4 = (1.0 - ibmx_h4 / (KmIbmxPde4 + ibmx_h4)) * PDE4_tot

    f_pde3_cav = r_pde34_frac * f_pde4_part * PDE4_tot / PDE3_tot
    f_pde3_cyt = 1.0 - f_pde3_cav

    c[50] = ibmx2 * f_pde2_cav * c[5]  # Concentration of PDE2 in Caveolar subspace
    c[51] = ibmx2 * f_pde2_eca * c[6]  # Concentration of PDE2 in Extracaveolar subspace
    c[52] = ibmx2 * f_pde2_cyt * c[7]  # Concentration of PDE2 in Cytosolic subspace

    c[53] = ibmx3 * f_pde3_cav * c[5]  # Concentration of PDE3 in Caveolar subspace
    c[54] = ibmx3 * f_pde3_cyt * c[7]  # Concentration of PDE3 in Cytosolic subspace

    c[55] = ibmx4 * f_pde4_cav * c[5]  # Concentration of PDE4 in Caveolar subspace
    c[56] = ibmx4 * f_pde4_eca * c[6]  # Concentration of PDE4 in Extracaveolar subspace
    c[57] = ibmx4 * f_pde4_cyt * c[7]  # Concentration of PDE4 in Cytosolic subspace

    # Adrenergic Receptor and G Protein Activation (pg 18 - 24 of Heijman supplementary document)
    beta_R_b1_tot = 0.85 * 0.025  # Total cellular beta-1 adrenergic receptor concentration
    beta_R_b2_tot = 0.15 * 0.025  # Total cellular beta-2 adrenergic receptor concentration
    c[58] = 224.0 * beta_R_b1_tot  # Total Gs protein concentration
    c[59] = 0.5  # Total Gi protein concentration

    c[60] = 0.5664  # Fraction of Gs proteins located in extracaveolar space
    c[61] = 0.0011071  # Fraction of Gs proteins located in caveolar subspace
    c[62] = 1.0 - c[61] - c[60]  # Fraction of Gs proteins located in cytoplasm

    c[63] = 0.85  # Fraction of Gi proteins located in caveolar subspace
    c[64] = 1.0 - c[63]  # Fraction of Gi proteins located in extracaveolar space

    f_Rb1_cav = 0.081161  # Fraction of beta-1 adrenergic receptors located in caveolar subspace
    f_Rb1_eca = 0.48744  # Fraction of beta-1 adrenergic receptors located in extra caveolar space
    f_Rb1_cyt = 1.0 - f_Rb1_cav - f_Rb1_eca  # Fraction of beta-1 adrenergic receptors located in cytoplasm

    f_Rb2_cav = 0.85  # Fraction of beta-2 adrenergic receptors located in caveolar subspace
    f_Rb2_eca = 1.0 - f_Rb2_cav  # Fraction of beta-2 adrenergic receptors located in extracaveolar space

    c[65] = 0.567  # beta-1 receptor / ligand low affinity constant
    c[66] = 2.449  # beta-1 receptor / G-protein affinity constant
    c[67] = 0.062  # beta-1 receptor / ligand high affinity constant

    c[68] = 1.053  # Phosph. beta-2 receptor / ligand low affinity constant
    c[69] = 0.012  # beta-2 receptor / ligand high affinity constant
    c[70] = 1.6655  # Phosph. beta-2 receptor / Gi-protein affinity constant
    c[71] = 1.8463  # beta-2 receptor / G-protein affinity constant
    c[72] = 1.053  # beta-2 receptor / ligand low affinity constant
    c[73] = 0.1  # Phosph. beta-2 receptor / ligand high affinity constant

    c[74] = 4.9054  # Activation rate for Gs by high affinity complex
    c[75] = 0.25945  # Activation rate for Gs by low affinity complex

    c[76] = 4.0  # Activation rate for Gi by high affinity complex
    c[77] = 0.05  # Activation rate for Gi by low affinity complex

    c[78] = 0.8  # Gs GTP to GDP hydrolysis constant
    c[79] = c[78]  # Gi GTP to GDP hydrolysis constant

    c[80] = 1210.0  # Reassociation rate for Gs subunits
    c[81] = c[80]  # Reassociation rate for Gi subunits

    rate_bds = 0.35

    c[82] = rate_bds * 0.0009833  # Rate for GRK dephosphorylation
    c[83] = rate_bds * 0.00133  # Rate for GRK dependent receptor desensitization

    c[84] = rate_bds * 0.0065  # Rate for (PKA dependent receptor) desensitization
    c[85] = 0.15629 * c[84]  # Rate for PKA dephosphorylation

    # beta_cav
    c[86] = 1.0
    c[87] = f_Rb1_cav * beta_R_b1_tot * c[5]  # Total concentration of beta1AR in the caveolar subspace
    c[88] = 1.0  # Ratio of beta2AR incorporated in RGs_Tot and LRGs_Tot
    c[89] = f_Rb2_cav * beta_R_b2_tot * c[5]  # Total concentration of beta2AR in the caveolar subspace
    c[90] = (c[73] + c[1]) * (c[68] + c[1]) / c[68]
    c[91] = c[71] * c[69] * c[65] * (c[67] + c[1]) * (c[72] + c[1])
    c[92] = c[65] * c[72] * (c[67] + c[1]) * (c[69] + c[1])
    c[93] = c[66] * c[71] * c[67] * c[69] * (c[65] + c[1]) * (c[72] + c[1])
    c[94] = c[66] * c[67] * c[72] * (c[69] + c[1]) * (c[65] + c[1])

    # beta_eca
    c[95] = 1.0
    c[96] = 1.0  # Ratio of beta2AR incorporated in RGs_Tot and LRGs_Tot
    c[97] = f_Rb2_eca * beta_R_b2_tot * c[6]  # Total concentration of beta2AR in the extracaveolar space
    c[98] = f_Rb1_eca * beta_R_b1_tot * c[6]  # Total concentration of beta1AR in the extracaveolar space
    c[99] = (c[73] + c[1]) * (c[68] + c[1]) / c[68]
    c[100] = c[66] * c[67] * c[72] * (c[69] + c[1]) * (c[65] + c[1])
    c[101] = c[66] * c[71] * c[67] * c[69] * (c[65] + c[1]) * (c[72] + c[1])
    c[102] = c[65] * c[72] * (c[67] + c[1]) * (c[69] + c[1])
    c[103] = c[71] * c[69] * c[65] * (c[67] + c[1]) * (c[72] + c[1])

    # beta_cyt
    c[104] = f_Rb1_cyt * beta_R_b1_tot * c[7]  # Total concentration of beta-1 AR in the cytoplasm
    c[105] = 1.0
    c[106] = (c[67] + c[1]) * (c[65] + c[1]) / c[65]

    # AC (see pg 24-26 of Heijman supplementary document)
    ATP = 5000.0  # Concentration of ATP
    KmATP = 315.0  # AC affinity for ATP
    AC_tot = 3.0 * beta_R_b1_tot  # Total cellular AC concentration

    f_AC47_eca = 0.16479  # Fraction of AC47 in extracaveolar space
    f_AC56_cav = 0.087459  # Fraction of AC56 located in caveolae
    f_AC56_AC47 = 1.0 / (1.0 + 0.35)  # Fraction of AC that is of type 5/6

    c[107] = 1.0043  # Hill coefficient for AC47 activation
    c[108] = 1.3574  # Hill coefficient for AC56 activation
    c[109] = 0.6623  # Hill coefficient for Gs/Gi interaction of AC56

    c[110] = 0.031544  # AC47 affinity for Gs
    c[111] = 0.0852  # AC56 affinity for Gs
    c[112] = 0.0465  # AC56 affinity for inhibition by Gi
    c[113] = 0.4824  # Gs-dependence of inactivation by Gi for AC56

    c[114] = 0.03135  # Basal AC47 activity
    c[115] = 0.037696  # Basal AC56 activity

    c[116] = 3.3757  # Amplification factor for AC47
    c[117] = 41.32  # Amplification factor for AC56

    c[118] = 0.8569  # Maximum reduction in Gi inhibition, by Gs

    c[119] = (1.0 - f_AC47_eca) * (1.0 - f_AC56_AC47) * AC_tot * c[7]  # Concentration of AC56 in the cytoplasm
    c[120] = f_AC56_cav * f_AC56_AC47 * AC_tot * c[5]  # Concentration of AC56 in the caveolar subspace

    c[121] = ATP / (KmATP + ATP)
    c[122] = f_AC47_eca * (1.0 - f_AC56_AC47) * AC_tot * c[6]  # Concentration of AC56 in the extracaveolar space
    c[123] = (1.0 - f_AC56_cav) * f_AC56_AC47 * AC_tot * c[7]  # Concentration of AC56 in the cytoplasm

    # Substrate Phosphorylation
    # iNaK
    c[124] = 0.015265  # Rate of INaK phosphorylation by PKA
    c[125] = 0.092455  # Rate of INaK dephosphorylation by phosphatases
    c[126] = 0.0011001  # Affinity of INaK for phosphorylation by PKA
    c[127] = 5.7392  # Affinity of INaK for dephosphorylation by phosphatases

    # iNa
    c[128] = 0.10988  # Affinity of INa for phosphorylation by PKA
    c[129] = 7.8605  # Affinity of INa for dephosphorylation by phosphatases
    c[130] = 0.01368  # Rate of INa phosphorylation by PKA
    c[131] = 0.052811  # Rate of INa dephosphorylation by phosphatases

    # iup
    c[132] = 0.11348  # Rate of PLB phosphorylation by PKA
    c[133] = 0.48302  # Rate of PLB dephosphorylation by phosphatases
    c[134] = 9.88539999999999992e-04  # Affinity of PLB for phosphorylation by PKA
    c[135] = 0.80737  # Affinity of PLB for dephosphorylation by phosphatases

    # iKur
    c[136] = 0.069537  # Rate of IKur phosphorylation by PKA
    c[137] = 0.317  # Rate of IKur dephosphorylation by phosphatases
    c[138] = 0.27623  # Affinity of IKur for phosphorylation by PKA
    c[139] = 0.002331  # Affinity of IKur for dephosphorylation by phosphatases

    # iKs
    c[140] = 0.16305  # Rate of IKs channel phosphorylation by PKA
    c[141] = 1.0542  # Rate of IKs channel dephosphorylation by phosphatases
    c[142] = 9.97940000000000003e-05  # Affinity of IKs channels for phosphorylation by PKA
    c[143] = 1.11470000000000002e-04  # Affinity of IKs channels for dephosphorylation by phosphatases

    M = 0.01  # Binding affinity between PKA and Yotiao
    iks_sig_L = 0.0001  # Binding affinity between PKA and Yotiao
    iks_sig_K = 0.01  # Binding affinity between PP1 and Yotiao
    Yotiao = 0.025  # Total concentration of Yotiao
    c[144] = 0.025  # Total concentration of IKs channels
    iks_sig_PKAf_sum = 1.0 + (Yotiao - c[11]) / M
    iks_sig_PKAf = (
        M / 2.0 * (√(iks_sig_PKAf_sum^2.0 + 4.0 * c[11] / M) - iks_sig_PKAf_sum)
    )  # Concentration of PKA not bound to AKAPs in the extracaveolar compartment
    iks_sig_PP1f_eca_sum = 1.0 + (Yotiao - c[35]) / iks_sig_K
    PP1f_eca = (
        iks_sig_K / 2.0 *
        (√(iks_sig_PP1f_eca_sum^2.0 + 4.0 * c[35] / iks_sig_K) - iks_sig_PP1f_eca_sum)
    )  # Concentration of PP1 not bound to AKAPs in the extracaveolar compartment
    iks_sig_IKsf_sum = 1.0 + (Yotiao - c[144]) / iks_sig_L
    IKsf = (
        iks_sig_L / 2.0 *
        (√(iks_sig_IKsf_sum^2.0 + 4.0 * c[144] / iks_sig_L) - iks_sig_IKsf_sum)
    )  # Concentration of IKs not bound to AKAPs in the extracaveolar compartment
    Yotiaof =
        (Yotiao - c[144] + IKsf) / ((1.0 + PP1f_eca / iks_sig_K) * (1.0 + iks_sig_PKAf / M))
    c[145] = (IKsf * Yotiaof * iks_sig_PKAf / (iks_sig_L * M))  # Concentration of IKs channels that have an AKAP and local PKA but no PP1 available
    c[146] = (c[145] * PP1f_eca / iks_sig_K)  # Concentration of IKs channels that have an AKAP, local PKA and PP1 available

    # Troponin
    c[147] = 0.10408  # Rate of Troponin phosphorylation by PKA
    c[148] = 0.052633  # Rate of Troponin dephosphorylation by phosphatases
    c[149] = 2.71430000000000008e-05  # Affinity of Troponin for phosphorylation by PKA
    c[150] = 0.26714  # Affinity of Troponin for dephosphorylation by phosphatases

    # RyR and iCaL  -  Substrates with A-Kinase Anchoring Protein (AKAP)
    c[151] = 0.125  # Total concentration of RyRs
    RyR_akap = 0.125  # Total concentration of RyR AKAP
    c[152] = 0.0025548  # Rate of RyR phosphorylation by PKA
    c[153] = 0.0038257  # Rate of RyR dephosphorylation by phosphatases
    c[154] = 6.62979999999999944e-05  # Affinity of RyR for phosphorylation by PKA
    c[155] = 0.043003  # Affinity of RyR for dephosphorylation by phosphatases
    Mr = 0.01  # Binding affinity between PKA and ICaL AKAP
    Lr = 0.0001  # Binding affinity between RyR and AKAP
    Kr = 0.01  # Binding affinity between PP1 and RyR AKAP

    akap_sig_RyRf_sum = 1.0 + (RyR_akap - c[151]) / Lr
    RyRf = (Lr / 2.0 * (√(akap_sig_RyRf_sum^2.0 + 4.0 * c[151] / Lr) - akap_sig_RyRf_sum))  # Caveolar concentration of free RyR

    c[156] = 0.025  # Total concentration of ICaL channels
    ICaL_akap = 0.025  # Total concentration of ICaL AKAP
    c[157] = 5.10090000000000044e-04  # Rate of ICaL channel phosphorylation by PKA
    c[158] = 0.0006903  # Rate of ICaL channel dephosphorylation by phosphatases
    c[159] = 1.27019999999999993e-06  # Affinity of ICaL channels for phosphorylation by PKA
    c[160] = 0.0063064  # Affinity of ICaL channels for dephosphorylation by phosphatases
    Mi = 0.01  # Binding affinity between PKA and ICaL AKAP
    Li = 0.0001  # Binding affinity between ICaL channel and AKAP
    Ki = 0.01  # Binding affinity between PP1 and ICaL AKAP
    akap_sig_ICaLf_sum = 1.0 + (ICaL_akap - c[156]) / Li
    ICaLf = (
        Li / 2.0 * (√(akap_sig_ICaLf_sum^2.0 + 4.0 * c[156] / Li) - akap_sig_ICaLf_sum)
    )  # Caveolar concentration of free ICaL

    akap_sig_PP1f_cav_b = ICaL_akap + RyR_akap + Ki + Kr - c[36]
    akap_sig_PP1f_cav_c = ICaL_akap * Kr + RyR_akap * Ki + Ki * Kr - c[36] * (Ki + Kr)
    akap_sig_PP1f_cav_d = c[36] * Ki * Kr
    akap_sig_PP1f_cav_rr = (
        -akap_sig_PP1f_cav_d / 27.0 * akap_sig_PP1f_cav_b^3.0 -
        akap_sig_PP1f_cav_b *
        akap_sig_PP1f_cav_b *
        akap_sig_PP1f_cav_c *
        akap_sig_PP1f_cav_c / 108.0 +
        akap_sig_PP1f_cav_b * akap_sig_PP1f_cav_c * akap_sig_PP1f_cav_d / 6.0 +
        akap_sig_PP1f_cav_c^3.0 / 27.0 +
        akap_sig_PP1f_cav_d * akap_sig_PP1f_cav_d / 4.0
    )
    akap_sig_PP1f_cav_yi = ((akap_sig_PP1f_cav_rr < 0.0) ? √(-akap_sig_PP1f_cav_rr) : 0.0)
    akap_sig_PP1f_cav_yr = (
        ((akap_sig_PP1f_cav_rr > 0.0) ? √(akap_sig_PP1f_cav_rr) : 0.0) +
        akap_sig_PP1f_cav_d / 2.0 +
        akap_sig_PP1f_cav_b * akap_sig_PP1f_cav_c / 6.0 - akap_sig_PP1f_cav_b^3.0 / 27.0
    )
    akap_sig_PP1f_cav_mag =
        (
            akap_sig_PP1f_cav_yr * akap_sig_PP1f_cav_yr +
            akap_sig_PP1f_cav_yi * akap_sig_PP1f_cav_yi
        )^(1.0 / 6.0)
    akap_sig_PP1f_cav_arg = atan(akap_sig_PP1f_cav_yi / akap_sig_PP1f_cav_yr) / 3.0
    akap_sig_PP1f_cav_x =
        (akap_sig_PP1f_cav_c / 3.0 - akap_sig_PP1f_cav_b * akap_sig_PP1f_cav_b / 9.0) /
        (akap_sig_PP1f_cav_mag * akap_sig_PP1f_cav_mag)

    PP1f_cav = (
        akap_sig_PP1f_cav_mag * cos(akap_sig_PP1f_cav_arg) * (1.0 - akap_sig_PP1f_cav_x) -
        akap_sig_PP1f_cav_b / 3.0
    )  # Caveolar concentration of free PP1
    akap_sig_PKAf_d = c[12] * Mi * Mr
    akap_sig_PKAf_b = ICaL_akap + RyR_akap + Mi + Mr - c[12]
    akap_sig_PKAf_c = ICaL_akap * Mr + RyR_akap * Mi + Mi * Mr - c[12] * (Mi + Mr)
    akap_sig_PKAf_rr = (
        -akap_sig_PKAf_d / 27.0 * akap_sig_PKAf_b^3.0 -
        akap_sig_PKAf_b * akap_sig_PKAf_b * akap_sig_PKAf_c * akap_sig_PKAf_c / 108.0 +
        akap_sig_PKAf_b * akap_sig_PKAf_c * akap_sig_PKAf_d / 6.0 +
        akap_sig_PKAf_c^3.0 / 27.0 +
        akap_sig_PKAf_d * akap_sig_PKAf_d / 4.0
    )
    akap_sig_PKAf_yr = (
        ((akap_sig_PKAf_rr > 0.0) ? √(akap_sig_PKAf_rr) : 0.0) +
        akap_sig_PKAf_d / 2.0 +
        akap_sig_PKAf_b * akap_sig_PKAf_c / 6.0 - akap_sig_PKAf_b^3.0 / 27.0
    )
    akap_sig_PKAf_yi = ((akap_sig_PKAf_rr < 0.0) ? √(-akap_sig_PKAf_rr) : 0.0)
    akap_sig_PKAf_mag =
        (
            akap_sig_PKAf_yr * akap_sig_PKAf_yr + akap_sig_PKAf_yi * akap_sig_PKAf_yi
        )^(1.0 / 6.0)
    akap_sig_PKAf_arg = atan(akap_sig_PKAf_yi / akap_sig_PKAf_yr) / 3.0
    akap_sig_PKAf_x =
        (akap_sig_PKAf_c / 3.0 - akap_sig_PKAf_b * akap_sig_PKAf_b / 9.0) /
        (akap_sig_PKAf_mag * akap_sig_PKAf_mag)

    akap_sig_PKAf = (
        akap_sig_PKAf_mag * cos(akap_sig_PKAf_arg) * (1.0 - akap_sig_PKAf_x) -
        akap_sig_PKAf_b / 3.0
    )  # Caveolar concentration of free PKA

    RyR_akapf =
        (RyR_akap - c[151] + RyRf) / ((PP1f_cav / Kr + 1.0) * (akap_sig_PKAf / Mr + 1.0))  # Caveolar concentration of free RyR AKAP
    c[161] = (RyRf * RyR_akapf * akap_sig_PKAf / (Lr * Mr))  # Concentration of RyR that have an AKAP, local PKA but no PP1 available
    c[162] = c[161] * PP1f_cav / Kr  # Concentration of RyR that have an AKAP, local PKA and PP1

    ICaL_akapf =
        (ICaL_akap - c[156] + ICaLf) / ((PP1f_cav / Ki + 1.0) * (akap_sig_PKAf / Mi + 1.0))  # Caveolar concentration of free ICaL AKAP
    c[163] = (ICaLf * ICaL_akapf * akap_sig_PKAf / (Li * Mi))  # Concentration of ICaL channels that have an AKAP, local PKA but no PP1 available
    c[164] = c[163] * PP1f_cav / Ki  # Concentration of ICaL channels that have an AKAP, local PKA and PP1

    # These are for electrophysiology.. check these to make sure
    # that the one for O'Hara is consisten
    # irel

    c[165] = 0.0329 + c[161] / c[151]

    # ical
    c[166] = 0.0269 + c[163] / c[156]

    # iks
    c[167] = 0.0306 + c[145] / c[144]

    # Return as named tuple for unpacking into @mtkmodel
    return (; (Symbol("c$i") => c[i] for i in 1:167)...)
end
