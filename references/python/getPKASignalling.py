import numpy as np


def get_pka_signalling(y: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Calculates the derivatives for the PKA signaling pathway model.

    This function is a Python conversion of the MATLAB code adapted from:
    https://github.com/JQXGong/Ohara-beta-adrenergic

    Reference:
    Quantitative analysis of variability in an integrated model of
    human ventricular electrophysiology and B-adrenergic signaling
    by Jingqi Q.X. Gong, Monica E. Susilo, Anna Sher, Cynthia J. Musante, Eric A. Sobie
    DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009

    Args:
        y (np.ndarray): A 1D NumPy array of the current state variables.
        c (np.ndarray): A 1D NumPy array of the model parameters/constants.

    Returns:
        np.ndarray: A 1D NumPy array containing the derivatives (ydot) of the state variables.
    """

    # Create derivatives vector
    ydot = np.zeros_like(y)

    # %% Unpack state variables from the y vector
    beta_cav_Gs_aGTP = y[0]  # 1: Gs_aGTP_CAV
    beta_eca_Gs_aGTP = y[1]  # 2: Gs_aGTP_ECAV
    beta_cyt_Gs_aGTP = y[2]  # 3: Gs_a_GTP_CYT
    beta_cav_Gs_bg = y[3]  # 4: Gs_bg_CAV
    beta_eca_Gs_bg = y[4]  # 5: Gs_bg_ECAV
    beta_cyt_Gs_bg = y[5]  # 6: Gs_bg_CYT
    beta_cav_Gs_aGDP = y[6]  # 7: Gs_aGDP_CAV
    beta_eca_Gs_aGDP = y[7]  # 8: Gs_aGDP_ECAV
    beta_cyt_Gs_aGDP = y[8]  # 9: Gs_aGDP_CYT

    cAMP_cav = y[9]  # 10: cAMP_CAVVV
    cAMP_eca = y[10]  # 11: cAMP_ECAV
    cAMP_cyt = y[11]  # 12: cAMP_CYT

    beta_cav_Rb1_pka_tot = y[12]  # 13: R_pkap_tot_CAV
    beta_eca_Rb1_pka_tot = y[13]  # 14: R_pkap_tot_ECAV
    beta_cyt_Rb1_pka_tot = y[14]  # 15: R_pkap_tot_CYT
    beta_cav_Rb1_grk_tot = y[15]  # 16: R_grkp_tot_CAV
    beta_eca_Rb1_grk_tot = y[16]  # 17: R_grkp_tot_ECAV
    beta_cyt_Rb1_grk_tot = y[17]  # 18: R_grkp_tot_CYT

    pka_cav_ARC = y[18]  # 19: RLC_CAV
    pka_cav_A2RC = y[19]  # 20: L2RC_CAV
    pka_cav_A2R = y[20]  # 21: L2R_CAV
    pka_cav_C = y[21]  # 22: C_CAV
    pka_cav_PKIC = y[22]  # 23: PKI_CAV
    pka_eca_ARC = y[23]  # 24: RLC_ECAV
    pka_eca_A2RC = y[24]  # 25: L2RC_ECAV
    pka_eca_A2R = y[25]  # 26: L2R_ECAV
    pka_eca_C = y[26]  # 27: C_ECAV
    pka_eca_PKIC = y[27]  # 28: PKI_ECAV
    pka_cyt_ARC = y[28]  # 29: RLC_CYT
    pka_cyt_A2RC = y[29]  # 30: L2RC_CYT
    pka_cyt_A2R = y[30]  # 31: L2R_CYT
    pka_cyt_C = y[31]  # 32: C_CYT
    pka_cyt_PKIC = y[32]  # 33: PKI_CYT

    PDE3_P_cav = y[33]  # 34: PDE3_P_CAV
    PDE3_P_cyt = y[34]  # 35: PDE3_P_CYT
    PDE4_P_cav = y[35]  # 36: PDE4_P_CAV
    PDE4_P_eca = y[36]  # 37: PDE4_P_ECAV
    PDE4_P_cyt = y[37]  # 38: PDE4_P_CYT

    inhib1_p = y[38]  # 39: Inhib1_P_CYT

    ICaLp = y[39]  # 40: fLCC_P
    IKsp = y[40]  # 41: fIKS_P
    iup_f_plb = y[41]  # 42: fPLB_P
    f_tni = y[42]  # 43: fTnI_P
    ina_f_ina = y[43]  # 44: fINa_P
    f_inak = y[44]  # 45: fINaK_P
    RyRp = y[45]  # 46: fRyR_P
    f_ikur = y[46]  # 47: fIKur_P

    # Clamp PKA fractions to be between [0,1] (using 0.0001 and 0.9999)
    # iup_f_plb
    if iup_f_plb < 0.0:
        iup_f_plb = 0.0001
    elif iup_f_plb > 1.0:
        iup_f_plb = 0.9999
    # f_tni
    if f_tni < 0.0:
        f_tni = 0.0001
    elif f_tni > 1.0:
        f_tni = 0.9999
    # ina_f_ina
    if ina_f_ina < 0.0:
        ina_f_ina = 0.0001
    elif ina_f_ina > 1.0:
        ina_f_ina = 0.9999
    # f_inak
    if f_inak < 0.0:
        f_inak = 0.0001
    elif f_inak > 1.0:
        f_inak = 0.9999
    # f_ikur
    if f_ikur < 0.0:
        f_ikur = 0.0001
    elif f_ikur > 1.0:
        f_ikur = 0.9999

    beta_cav_Rb2_pka_tot = y[47]  # 48: Rb2_pkap_tot_CAV
    beta_cav_Rb2_grk_tot = y[48]  # 49: Rb2_grkp_tot_CAV
    beta_cav_Gi_aGTP = y[49]  # 50: Gi_aGTP_CAV
    beta_cav_Gi_bg = y[50]  # 51: Gi_bg_CAV
    beta_cav_Gi_aGDP = y[51]  # 52: Gi_aGDP_CAV
    beta_eca_Rb2_pka_tot = y[52]  # 53: Rb2_pkap_tot_ECAV
    beta_eca_Rb2_grk_tot = y[53]  # 54: Rb2_grkp_tot_ECAV
    beta_eca_Gi_aGTP = y[54]  # 55: Gi_aGTP_ECAV
    beta_eca_Gi_bg = y[55]  # 56: Gi_bg_ECAV
    beta_eca_Gi_aGDP = y[56]  # 57: Gi_aGDP_ECAV

    # %% CAVEOLAR COMPARTMENT %%
    # --------------------------------------------------------------------------
    beta_cav_Rb1_np_tot = c[86] - beta_cav_Rb1_pka_tot - beta_cav_Rb1_grk_tot
    beta_cav_Rb2_np_tot = c[88] - beta_cav_Rb2_pka_tot - beta_cav_Rb2_grk_tot
    beta_cav_Gi_abg = c[62] * c[58] * c[4] - beta_cav_Gi_aGTP - beta_cav_Gi_aGDP
    beta_cav_Gs_abg = c[60] * c[57] * c[4] - beta_cav_Gs_aGTP - beta_cav_Gs_aGDP

    beta_cav_Gs_f_d = beta_cav_Gs_abg * c[92] / c[91]
    beta_cav_Gs_f_b = (
        (c[93] + c[90]) / c[91]
        + beta_cav_Rb1_np_tot
        + beta_cav_Rb2_np_tot
        - beta_cav_Gs_abg
    )
    beta_cav_Gs_f_c = (
        c[90] * (beta_cav_Rb1_np_tot - beta_cav_Gs_abg)
        + c[93] * (beta_cav_Rb2_np_tot - beta_cav_Gs_abg)
        + c[92]
    ) / c[91]
    beta_cav_Gs_f_rr = (
        -beta_cav_Gs_f_d / 27.0 * beta_cav_Gs_f_b**3.0
        - beta_cav_Gs_f_b * beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_c / 108.0
        + beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_d / 6.0
        + beta_cav_Gs_f_c**3.0 / 27.0
        + beta_cav_Gs_f_d * beta_cav_Gs_f_d / 4.0
    )
    beta_cav_Gs_f_yr = (
        (np.sqrt(beta_cav_Gs_f_rr) if (beta_cav_Gs_f_rr > 0.0) else 0.0)
        + beta_cav_Gs_f_d / 2.0
        + beta_cav_Gs_f_b * beta_cav_Gs_f_c / 6.0
        - beta_cav_Gs_f_b**3.0 / 27.0
    )
    beta_cav_Gs_f_yi = np.sqrt(-beta_cav_Gs_f_rr) if (beta_cav_Gs_f_rr < 0.0) else 0.0
    beta_cav_Gs_f_mag = (
        beta_cav_Gs_f_yr * beta_cav_Gs_f_yr + beta_cav_Gs_f_yi * beta_cav_Gs_f_yi
    ) ** (1.0 / 6.0)
    beta_cav_Gs_f_arg = (
        np.arctan(beta_cav_Gs_f_yi / beta_cav_Gs_f_yr) / 3.0
        if beta_cav_Gs_f_yr != 0
        else 0
    )
    beta_cav_Gs_f_x = (
        beta_cav_Gs_f_c / 3.0 - beta_cav_Gs_f_b * beta_cav_Gs_f_b / 9.0
    ) / (beta_cav_Gs_f_mag * beta_cav_Gs_f_mag)
    beta_cav_Gs_f_r = (
        beta_cav_Gs_f_mag * np.cos(beta_cav_Gs_f_arg) * (1.0 - beta_cav_Gs_f_x)
        - beta_cav_Gs_f_b / 3.0
    )
    beta_cav_Gs_f_i = (
        beta_cav_Gs_f_mag * np.sin(beta_cav_Gs_f_arg) * (1.0 + beta_cav_Gs_f_x)
    )

    beta_cav_Gs_f = np.sqrt(
        beta_cav_Gs_f_r * beta_cav_Gs_f_r + beta_cav_Gs_f_i * beta_cav_Gs_f_i
    )
    beta_cav_Rb1_f = beta_cav_Rb1_np_tot / (
        1.0 + c[0] / c[64] + beta_cav_Gs_f * (c[66] + c[0]) / (c[65] * c[66])
    )
    beta_cav_LRb1 = c[0] * beta_cav_Rb1_f / c[64]
    beta_cav_LRb1Gs = c[0] * beta_cav_Rb1_f * beta_cav_Gs_f / (c[65] * c[66])
    beta_cav_Rb2_f = beta_cav_Rb2_np_tot / (
        1.0 + c[0] / c[71] + beta_cav_Gs_f * (c[68] + c[0]) / (c[70] * c[68])
    )
    beta_cav_LRb2 = c[0] * beta_cav_Rb2_f / c[71]
    beta_cav_Rb2Gs = beta_cav_Rb2_f * beta_cav_Gs_f / c[70]
    beta_cav_Rb1Gs = beta_cav_Rb1_f * beta_cav_Gs_f / c[65]
    beta_cav_LRb2Gs = c[0] * beta_cav_Rb2_f * beta_cav_Gs_f / (c[70] * c[68])

    ydot[12] = 0.001 * (
        c[83] * pka_cav_C * beta_cav_Rb1_np_tot - c[84] * beta_cav_Rb1_pka_tot
    )
    ydot[15] = 0.001 * (
        c[82] * c[85] * (beta_cav_LRb1 + beta_cav_LRb1Gs) - c[81] * beta_cav_Rb1_grk_tot
    )
    ydot[47] = 0.001 * (
        c[83] * pka_cav_C * beta_cav_Rb2_np_tot - c[84] * beta_cav_Rb2_pka_tot
    )
    ydot[48] = 0.001 * (
        c[82] * c[85] * (beta_cav_LRb2 + beta_cav_LRb2Gs) - c[81] * beta_cav_Rb2_grk_tot
    )

    # %% EXTRACAVEOLAR COMPARTMENT %%
    # --------------------------------------------------------------------------
    beta_eca_Gs_abg = c[59] * c[57] * c[5] - beta_eca_Gs_aGTP - beta_eca_Gs_aGDP
    beta_eca_Rb2_np_tot = c[96] - beta_eca_Rb2_pka_tot - beta_eca_Rb2_grk_tot
    beta_eca_Rb1_np_tot = c[97] - beta_eca_Rb1_pka_tot - beta_eca_Rb1_grk_tot
    beta_eca_Gs_f_d = beta_eca_Gs_abg * c[100] / c[101]
    beta_eca_Gs_f_c = (
        c[102] * (beta_eca_Rb1_np_tot - beta_eca_Gs_abg)
        + c[99] * (beta_eca_Rb2_np_tot - beta_eca_Gs_abg)
        + c[100]
    ) / c[101]
    beta_eca_Gs_f_b = (
        (c[99] + c[102]) / c[101]
        + beta_eca_Rb1_np_tot
        + beta_eca_Rb2_np_tot
        - beta_eca_Gs_abg
    )
    beta_eca_Gs_f_rr = (
        -beta_eca_Gs_f_d / 27.0 * beta_eca_Gs_f_b**3.0
        - beta_eca_Gs_f_b * beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_c / 108.0
        + beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_d / 6.0
        + beta_eca_Gs_f_c**3.0 / 27.0
        + beta_eca_Gs_f_d * beta_eca_Gs_f_d / 4.0
    )
    beta_eca_Gs_f_yi = np.sqrt(-beta_eca_Gs_f_rr) if (beta_eca_Gs_f_rr < 0.0) else 0.0
    beta_eca_Gs_f_yr = (
        (np.sqrt(beta_eca_Gs_f_rr) if (beta_eca_Gs_f_rr > 0.0) else 0.0)
        + beta_eca_Gs_f_d / 2.0
        + beta_eca_Gs_f_b * beta_eca_Gs_f_c / 6.0
        - beta_eca_Gs_f_b**3.0 / 27.0
    )
    beta_eca_Gs_f_mag = (
        beta_eca_Gs_f_yr * beta_eca_Gs_f_yr + beta_eca_Gs_f_yi * beta_eca_Gs_f_yi
    ) ** (1.0 / 6.0)
    beta_eca_Gs_f_arg = (
        np.arctan(beta_eca_Gs_f_yi / beta_eca_Gs_f_yr) / 3.0
        if beta_eca_Gs_f_yr != 0
        else 0
    )
    beta_eca_Gs_f_x = (
        beta_eca_Gs_f_c / 3.0 - beta_eca_Gs_f_b * beta_eca_Gs_f_b / 9.0
    ) / (beta_eca_Gs_f_mag * beta_eca_Gs_f_mag)
    beta_eca_Gs_f_i = (
        beta_eca_Gs_f_mag * np.sin(beta_eca_Gs_f_arg) * (1.0 + beta_eca_Gs_f_x)
    )
    beta_eca_Gs_f_r = (
        beta_eca_Gs_f_mag * np.cos(beta_eca_Gs_f_arg) * (1.0 - beta_eca_Gs_f_x)
        - beta_eca_Gs_f_b / 3.0
    )
    beta_eca_Gs_f = np.sqrt(
        beta_eca_Gs_f_r * beta_eca_Gs_f_r + beta_eca_Gs_f_i * beta_eca_Gs_f_i
    )
    beta_eca_Rb1_f = beta_eca_Rb1_np_tot / (
        1.0 + c[0] / c[64] + beta_eca_Gs_f * (c[66] + c[0]) / (c[65] * c[66])
    )
    beta_eca_Rb2_f = beta_eca_Rb2_np_tot / (
        1.0 + c[0] / c[71] + beta_eca_Gs_f * (c[68] + c[0]) / (c[70] * c[68])
    )
    beta_eca_LRb1 = c[0] * beta_eca_Rb1_f / c[64]
    beta_eca_LRb2 = c[0] * beta_eca_Rb2_f / c[71]
    beta_eca_LRb2Gs = c[0] * beta_eca_Rb2_f * beta_eca_Gs_f / (c[70] * c[68])
    beta_eca_LRb1Gs = c[0] * beta_eca_Rb1_f * beta_eca_Gs_f / (c[65] * c[66])
    beta_eca_Rb2Gs = beta_eca_Rb2_f * beta_eca_Gs_f / c[70]
    beta_eca_Rb1Gs = beta_eca_Rb1_f * beta_eca_Gs_f / c[65]
    beta_eca_RGs_tot = beta_eca_Rb1Gs + c[95] * beta_eca_Rb2Gs
    beta_eca_LRGs_tot = beta_eca_LRb1Gs + c[95] * beta_eca_LRb2Gs
    beta_eca_Gi_abg = c[63] * c[58] * c[5] - beta_eca_Gi_aGTP - beta_eca_Gi_aGDP

    beta_eca_Rb2_pka_f_c = -beta_eca_Rb2_pka_tot * c[69] * c[72]
    beta_eca_Rb2_pka_f_b = (
        beta_eca_Gi_abg * (c[0] + c[72])
        - beta_eca_Rb2_pka_tot * (c[72] + c[0])
        + c[69] * c[72] * (1.0 + c[0] / c[67])
    )
    sqrt_term = (
        beta_eca_Rb2_pka_f_b * beta_eca_Rb2_pka_f_b - 4.0 * c[98] * beta_eca_Rb2_pka_f_c
    )
    beta_eca_Rb2_pka_f = (
        -beta_eca_Rb2_pka_f_b + np.sqrt(sqrt_term if sqrt_term > 0 else 0)
    ) / (2.0 * c[98])

    ydot[13] = 0.001 * (
        c[83] * pka_eca_C * beta_eca_Rb1_np_tot - c[84] * beta_eca_Rb1_pka_tot
    )
    ydot[16] = 0.001 * (
        c[82] * c[94] * (beta_eca_LRb1 + beta_eca_LRb1Gs) - c[81] * beta_eca_Rb1_grk_tot
    )
    ydot[52] = 0.001 * (
        c[83] * pka_eca_C * beta_eca_Rb2_np_tot - c[84] * beta_eca_Rb2_pka_tot
    )
    ydot[53] = 0.001 * (
        c[82] * c[94] * (beta_eca_LRb2 + beta_eca_LRb2Gs) - c[81] * beta_eca_Rb2_grk_tot
    )

    # %% CYTOPLASM %%
    # --------------------------------------------------------------------------
    beta_cyt_Gs_abg = c[61] * c[57] * c[6] - beta_cyt_Gs_aGTP - beta_cyt_Gs_aGDP
    beta_cyt_Rb1_np_tot = c[103] - beta_cyt_Rb1_pka_tot - beta_cyt_Rb1_grk_tot

    beta_cyt_Rb1_np_f_b = (
        beta_cyt_Gs_abg * (c[66] + c[0])
        - beta_cyt_Rb1_np_tot * (c[66] + c[0])
        + c[65] * c[66] * (1.0 + c[0] / c[64])
    )
    beta_cyt_Rb1_np_f_c = -beta_cyt_Rb1_np_tot * c[66] * c[65]
    sqrt_term = (
        beta_cyt_Rb1_np_f_b * beta_cyt_Rb1_np_f_b - 4.0 * c[105] * beta_cyt_Rb1_np_f_c
    )
    Rb1_np_f = (-beta_cyt_Rb1_np_f_b + np.sqrt(sqrt_term if sqrt_term > 0 else 0)) / (
        2.0 * c[105]
    )
    beta_cyt_Gs_f = beta_cyt_Gs_abg / (1.0 + Rb1_np_f / c[65] * (1.0 + c[0] / c[66]))
    LRb1_np = c[0] * Rb1_np_f / c[64]
    LRb1Gs_np = c[0] * Rb1_np_f * beta_cyt_Gs_f / (c[65] * c[66])
    Rb1Gs_np = beta_cyt_Gs_f * Rb1_np_f / c[65]

    ydot[14] = 0.001 * (
        c[83] * pka_cyt_C * beta_cyt_Rb1_np_tot - c[84] * beta_cyt_Rb1_pka_tot
    )
    ydot[17] = 0.001 * (
        c[82] * c[104] * (LRb1_np + LRb1Gs_np) - c[81] * beta_cyt_Rb1_grk_tot
    )

    # %% G-Protein Activation %%
    # --------------------------------------------------------------------------
    beta_cav_RGs_tot = beta_cav_Rb1Gs + c[87] * beta_cav_Rb2Gs
    beta_cav_LRGs_tot = beta_cav_LRb1Gs + c[87] * beta_cav_LRb2Gs
    beta_cav_Rb2_pka_f_c = -beta_cav_Rb2_pka_tot * c[69] * c[72]
    beta_cav_Rb2_pka_f_b = (
        beta_cav_Gi_abg * (c[0] + c[72])
        - beta_cav_Rb2_pka_tot * (c[72] + c[0])
        + c[69] * c[72] * (1.0 + c[0] / c[67])
    )
    sqrt_term = (
        beta_cav_Rb2_pka_f_b * beta_cav_Rb2_pka_f_b - 4.0 * c[89] * beta_cav_Rb2_pka_f_c
    )
    beta_cav_Rb2_pka_f = (
        -beta_cav_Rb2_pka_f_b + np.sqrt(sqrt_term if sqrt_term > 0 else 0)
    ) / (2.0 * c[89])
    beta_cav_Gi_f = beta_cav_Gi_abg / (
        1.0 + beta_cav_Rb2_pka_f / c[69] * (1.0 + c[0] / c[72])
    )
    beta_cav_Rb2Gi = beta_cav_Rb2_pka_f * beta_cav_Gi_f / c[69]
    beta_cav_LRb2Gi = beta_cav_Rb2Gi * c[0] / c[72]
    beta_eca_Gi_f = beta_eca_Gi_abg / (
        1.0 + beta_eca_Rb2_pka_f / c[69] * (1.0 + c[0] / c[72])
    )
    beta_eca_Rb2Gi = beta_eca_Rb2_pka_f * beta_eca_Gi_f / c[69]
    beta_eca_LRb2Gi = c[0] / c[72] * beta_eca_Rb2Gi

    ydot[0] = 0.001 * (
        c[74] * beta_cav_RGs_tot + c[73] * beta_cav_LRGs_tot - c[77] * beta_cav_Gs_aGTP
    )
    ydot[49] = 0.001 * (
        c[76] * beta_cav_Rb2Gi + c[75] * beta_cav_LRb2Gi - c[78] * beta_cav_Gi_aGTP
    )
    ydot[1] = 0.001 * (
        c[74] * beta_eca_RGs_tot + c[73] * beta_eca_LRGs_tot - c[77] * beta_eca_Gs_aGTP
    )
    ydot[54] = 0.001 * (
        c[76] * beta_eca_Rb2Gi + c[75] * beta_eca_LRb2Gi - c[78] * beta_eca_Gi_aGTP
    )
    ydot[2] = 0.001 * (c[74] * Rb1Gs_np + c[73] * LRb1Gs_np - c[77] * beta_cyt_Gs_aGTP)

    ydot[3] = 0.001 * (
        c[74] * beta_cav_RGs_tot
        + c[73] * beta_cav_LRGs_tot
        - c[79] * beta_cav_Gs_bg * beta_cav_Gs_aGDP
    )
    ydot[50] = 0.001 * (
        c[76] * beta_cav_Rb2Gi
        + c[75] * beta_cav_LRb2Gi
        - c[80] * beta_cav_Gi_bg * beta_cav_Gi_aGDP
    )
    ydot[4] = 0.001 * (
        c[74] * beta_eca_RGs_tot
        + c[73] * beta_eca_LRGs_tot
        - c[79] * beta_eca_Gs_bg * beta_eca_Gs_aGDP
    )
    ydot[55] = 0.001 * (
        c[76] * beta_eca_Rb2Gi
        + c[75] * beta_eca_LRb2Gi
        - c[80] * beta_eca_Gi_bg * beta_eca_Gi_aGDP
    )
    ydot[5] = 0.001 * (
        c[74] * Rb1Gs_np + c[73] * LRb1Gs_np - c[79] * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP
    )

    ydot[6] = 0.001 * (
        c[77] * beta_cav_Gs_aGTP - c[79] * beta_cav_Gs_bg * beta_cav_Gs_aGDP
    )
    ydot[51] = 0.001 * (
        c[78] * beta_cav_Gi_aGTP - c[80] * beta_cav_Gi_bg * beta_cav_Gi_aGDP
    )
    ydot[7] = 0.001 * (
        c[77] * beta_eca_Gs_aGTP - c[79] * beta_eca_Gs_bg * beta_eca_Gs_aGDP
    )
    ydot[56] = 0.001 * (
        c[78] * beta_eca_Gi_aGTP - c[80] * beta_eca_Gi_bg * beta_eca_Gi_aGDP
    )
    ydot[8] = 0.001 * (
        c[77] * beta_cyt_Gs_aGTP - c[79] * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP
    )

    # %% PKA Activation %%
    # --------------------------------------------------------------------------
    pka_cav_RCf = c[11] - pka_cav_ARC - pka_cav_A2RC - pka_cav_A2R

    ydot[18] = 0.001 * (
        c[18] * pka_cav_RCf * cAMP_cav
        - c[21] * pka_cav_ARC
        - c[19] * pka_cav_ARC * cAMP_cav
        + c[22] * pka_cav_A2RC
    )
    ydot[19] = 0.001 * (
        c[19] * pka_cav_ARC * cAMP_cav
        - (c[22] + c[20]) * pka_cav_A2RC
        + c[23] * pka_cav_A2R * pka_cav_C
    )
    ydot[20] = 0.001 * (c[20] * pka_cav_A2RC - c[23] * pka_cav_A2R * pka_cav_C)
    ydot[21] = 0.001 * (
        c[20] * pka_cav_A2RC
        - c[23] * pka_cav_A2R * pka_cav_C
        + c[17] * pka_cav_PKIC
        - c[16] * (c[13] - pka_cav_PKIC) * pka_cav_C
    )
    ydot[22] = 0.001 * (
        c[16] * (c[13] - pka_cav_PKIC) * pka_cav_C - c[17] * pka_cav_PKIC
    )

    pka_eca_RCf = c[10] - pka_eca_ARC - pka_eca_A2RC - pka_eca_A2R

    ydot[23] = 0.001 * (
        c[18] * pka_eca_RCf * cAMP_eca
        - c[24] * pka_eca_ARC
        - c[19] * pka_eca_ARC * cAMP_eca
        + c[25] * pka_eca_A2RC
    )
    ydot[24] = 0.001 * (
        c[19] * pka_eca_ARC * cAMP_eca
        - (c[25] + c[20]) * pka_eca_A2RC
        + c[26] * pka_eca_A2R * pka_eca_C
    )
    ydot[25] = 0.001 * (c[20] * pka_eca_A2RC - c[26] * pka_eca_A2R * pka_eca_C)
    ydot[26] = 0.001 * (
        c[20] * pka_eca_A2RC
        - c[26] * pka_eca_A2R * pka_eca_C
        + c[17] * pka_eca_PKIC
        - c[16] * (c[14] - pka_eca_PKIC) * pka_eca_C
    )
    ydot[27] = 0.001 * (
        c[16] * (c[14] - pka_eca_PKIC) * pka_eca_C - c[17] * pka_eca_PKIC
    )

    pka_cyt_RCf = c[12] - pka_cyt_ARC - pka_cyt_A2RC - pka_cyt_A2R

    ydot[28] = 0.001 * (
        c[18] * pka_cyt_RCf * cAMP_cyt
        - c[27] * pka_cyt_ARC
        - c[19] * pka_cyt_ARC * cAMP_cyt
        + c[28] * pka_cyt_A2RC
    )
    ydot[29] = 0.001 * (
        c[19] * pka_cyt_ARC * cAMP_cyt
        - (c[28] + c[20]) * pka_cyt_A2RC
        + c[29] * pka_cyt_A2R * pka_cyt_C
    )
    ydot[30] = 0.001 * (c[20] * pka_cyt_A2RC - c[29] * pka_cyt_A2R * pka_cyt_C)
    ydot[31] = 0.001 * (
        c[20] * pka_cyt_A2RC
        - c[29] * pka_cyt_A2R * pka_cyt_C
        + c[17] * pka_cyt_PKIC
        - c[16] * (c[15] - pka_cyt_PKIC) * pka_cyt_C
    )
    ydot[32] = 0.001 * (
        c[16] * (c[15] - pka_cyt_PKIC) * pka_cyt_C - c[17] * pka_cyt_PKIC
    )

    # %% cAMP Dynamics %%
    # --------------------------------------------------------------------------
    pka_cav_dcAMP = (
        -c[18] * pka_cav_RCf * cAMP_cav
        + c[21] * pka_cav_ARC
        - c[19] * pka_cav_ARC * cAMP_cav
        + c[22] * pka_cav_A2RC
    )
    pka_eca_dcAMP = (
        -c[18] * pka_eca_RCf * cAMP_eca
        + c[24] * pka_eca_ARC
        - c[19] * pka_eca_ARC * cAMP_eca
        + c[25] * pka_eca_A2RC
    )
    pka_cyt_dcAMP = (
        -c[18] * pka_cyt_RCf * cAMP_cyt
        + c[27] * pka_cyt_ARC
        - c[19] * pka_cyt_ARC * cAMP_cyt
        + c[28] * pka_cyt_A2RC
    )

    dcAMP_PDE2_cyt = c[51] * c[40] / (1.0 + c[43] / cAMP_cyt)
    dcAMP_PDE2_eca = c[50] * c[40] / (1.0 + c[43] / cAMP_eca)
    dcAMP_PDE2_cav = c[49] * c[40] / (1.0 + c[43] / cAMP_cav)
    dcAMP_PDE3_cav = (
        (c[52] + (c[46] - 1.0) * PDE3_P_cav) * c[41] / (1.0 + c[44] / cAMP_cav)
    )
    dcAMP_PDE4_cyt = (
        (c[56] + (c[46] - 1.0) * PDE4_P_cyt) * c[42] / (1.0 + c[45] / cAMP_cyt)
    )
    dcAMP_PDE4_eca = (
        (c[55] + (c[46] - 1.0) * PDE4_P_eca) * c[42] / (1.0 + c[45] / cAMP_eca)
    )
    dcAMP_PDE4_cav = (
        (c[54] + (c[46] - 1.0) * PDE4_P_cav) * c[42] / (1.0 + c[45] / cAMP_cav)
    )
    dcAMP_PDE3_cyt = (
        (c[53] + (c[46] - 1.0) * PDE3_P_cyt) * c[41] / (1.0 + c[44] / cAMP_cyt)
    )

    camp_cAMP_cyt_pde = dcAMP_PDE2_cyt + dcAMP_PDE3_cyt + dcAMP_PDE4_cyt
    camp_cAMP_cyt_j1 = c[8] * (cAMP_cav - cAMP_cyt) / c[3]
    camp_cAMP_cyt_j2 = c[9] * (cAMP_eca - cAMP_cyt) / c[3]

    camp_cAMP_eca_pde = dcAMP_PDE2_eca + dcAMP_PDE4_eca
    camp_cAMP_eca_j2 = c[9] * (cAMP_eca - cAMP_cyt) / c[2]
    camp_cAMP_eca_j1 = c[7] * (cAMP_cav - cAMP_eca) / c[2]

    camp_cAMP_cav_pde = dcAMP_PDE2_cav + dcAMP_PDE3_cav + dcAMP_PDE4_cav
    camp_cAMP_cav_j2 = c[8] * (cAMP_cav - cAMP_cyt) / c[1]
    camp_cAMP_cav_j1 = c[7] * (cAMP_cav - cAMP_eca) / c[1]

    ac_kAC47_cyt_gsa = beta_cyt_Gs_aGTP ** c[106]
    kAC47_cyt = c[115] * (c[113] + ac_kAC47_cyt_gsa / (c[109] + ac_kAC47_cyt_gsa))
    ac_kAC56_cav_gsa = beta_cav_Gs_aGTP ** c[107]
    gsi = beta_cav_Gs_aGTP ** c[108]
    kAC56_cav = (
        c[116]
        * (c[114] + ac_kAC56_cav_gsa / (c[110] + ac_kAC56_cav_gsa))
        * (
            1.0
            - (1.0 - c[117] * gsi / (c[112] + gsi))
            * beta_cav_Gi_bg
            / (c[111] + beta_cav_Gi_bg)
        )
    )
    ac_kAC47_eca_gsa = beta_eca_Gs_aGTP ** c[106]
    kAC47_eca = c[115] * (c[113] + ac_kAC47_eca_gsa / (c[109] + ac_kAC47_eca_gsa))
    ac_kAC56_cyt_gsa = beta_cyt_Gs_aGTP ** c[107]
    kAC56_cyt = c[116] * (c[114] + ac_kAC56_cyt_gsa / (c[110] + ac_kAC56_cyt_gsa))

    dcAMP_AC47_cyt = kAC47_cyt * c[118] * c[120]
    dcAMP_AC56_cyt = kAC56_cyt * c[122] * c[120]
    dcAMP_AC56_cav = kAC56_cav * c[119] * c[120]
    dcAMP_AC47_eca = kAC47_eca * c[121] * c[120]

    ydot[9] = 0.001 * (
        pka_cav_dcAMP
        + dcAMP_AC56_cav
        - camp_cAMP_cav_pde
        - camp_cAMP_cav_j1
        - camp_cAMP_cav_j2
    )
    ydot[10] = 0.001 * (
        pka_eca_dcAMP
        + dcAMP_AC47_eca
        - camp_cAMP_eca_pde
        + camp_cAMP_eca_j1
        - camp_cAMP_eca_j2
    )
    ydot[11] = 0.001 * (
        pka_cyt_dcAMP
        + dcAMP_AC47_cyt
        + dcAMP_AC56_cyt
        - camp_cAMP_cyt_pde
        + camp_cAMP_cyt_j1
        + camp_cAMP_cyt_j2
    )

    # %% PDE Phosphorylation %%
    # --------------------------------------------------------------------------
    ydot[33] = 0.001 * (c[47] * pka_cav_C * (c[52] - PDE3_P_cav) - c[48] * PDE3_P_cav)
    ydot[34] = 0.001 * (c[47] * pka_cyt_C * (c[53] - PDE3_P_cyt) - c[48] * PDE3_P_cyt)
    ydot[35] = 0.001 * (c[47] * pka_cav_C * (c[54] - PDE4_P_cav) - c[48] * PDE4_P_cav)
    ydot[36] = 0.001 * (c[47] * pka_eca_C * (c[55] - PDE4_P_eca) - c[48] * PDE4_P_eca)
    ydot[37] = 0.001 * (c[47] * pka_cyt_C * (c[56] - PDE4_P_cyt) - c[48] * PDE4_P_cyt)

    # %% PP1 Inhibition %%
    # --------------------------------------------------------------------------
    pp1_PP1f_cyt_sum = c[37] - c[36] + inhib1_p
    PP1f_cyt = 0.5 * (
        np.sqrt(pp1_PP1f_cyt_sum**2.0 + 4.0 * c[37] * c[36]) - pp1_PP1f_cyt_sum
    )
    di = c[38] - inhib1_p
    ydot[38] = 0.001 * (
        c[30] * pka_cyt_C * di / (c[32] + di)
        - c[31] * c[39] * inhib1_p / (c[33] + inhib1_p)
    )

    # %% Channel Phosphorylation %%
    # --------------------------------------------------------------------------
    # Substrates without AKAP
    ydot[41] = 0.001 * (
        c[131] * pka_cyt_C * (1.0 - iup_f_plb) / (c[133] + 1.0 - iup_f_plb)
        - c[132] * PP1f_cyt * iup_f_plb / (c[134] + iup_f_plb)
    )
    ydot[42] = 0.001 * (
        c[146] * pka_cyt_C * (1.0 - f_tni) / (c[148] + 1.0 - f_tni)
        - c[147] * c[39] * f_tni / (c[149] + f_tni)
    )
    ydot[43] = 0.001 * (
        c[129] * pka_cav_C * (1.0 - ina_f_ina) / (c[127] + 1.0 - ina_f_ina)
        - c[130] * c[35] * ina_f_ina / (c[128] + ina_f_ina)
    )
    ydot[44] = 0.001 * (
        c[123] * pka_cav_C * (1.0 - f_inak) / (c[125] + 1.0 - f_inak)
        - c[124] * c[35] * f_inak / (c[126] + f_inak)
    )
    ydot[46] = 0.001 * (
        c[135] * pka_eca_C * (1.0 - f_ikur) / (c[137] + 1.0 - f_ikur)
        - c[136] * c[34] * f_ikur / (c[138] + f_ikur)
    )

    # Substrates with AKAP
    iks_sig_IKsp_dif = c[145] - IKsp
    ydot[40] = 0.001 * (
        c[139] * pka_eca_C * iks_sig_IKsp_dif / (c[141] + iks_sig_IKsp_dif)
        - c[140] * c[34] * IKsp / (c[142] + IKsp)
    )
    akap_sig_RyRp_dif = c[161] - RyRp
    ydot[45] = 0.001 * (
        c[151] * pka_cav_C * akap_sig_RyRp_dif / (c[153] + akap_sig_RyRp_dif)
        - c[152] * c[35] * RyRp / (c[154] + RyRp)
    )
    akap_sig_ICaLp_dif = c[163] - ICaLp
    ydot[39] = 0.001 * (
        c[156] * pka_cav_C * akap_sig_ICaLp_dif / (c[158] + akap_sig_ICaLp_dif)
        - c[157] * c[35] * ICaLp / (c[159] + ICaLp)
    )

    return ydot
