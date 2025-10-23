import numpy as np


def get_effective_fraction(
    y: np.ndarray, c: np.ndarray, output: np.ndarray, original_output=False
) -> np.ndarray:
    """
    Calculates the effective fraction of phosphorylated substrates.

    This function is a Python conversion of the MATLAB code adapted from:
    https://github.com/JQXGong/Ohara-beta-adrenergic

    Reference:
    Quantitative analysis of variability in an integrated model of
    human ventricular electrophysiology and B-adrenergic signaling
    by Jingqi Q.X. Gong, Monica E. Susilo, Anna Sher, Cynthia J. Musante, Eric A. Sobie
    DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009

    Args:
        y (list or np.ndarray): A state vector containing model variables.
        c (list or np.ndarray): A vector of model parameters.

    Returns:
        tuple: A tuple containing the effective fractions for ICaL, IKs, PLB, TnI,
               INa, INaK, RyR, and IKur.
    """

    # -----------------------------------------------------------------
    # ICaL
    # -----------------------------------------------------------------
    icalp = y[39]
    # Fraction of phosphorylated ICaL channels
    fp_ical = (icalp + c[162]) / c[155]
    fp_ical = max(0.0001, fp_ical)
    fp_ical = min(0.9999, fp_ical)

    # Effective fraction of phosphorylated ICaL channels
    ical_f_hat_val = (fp_ical - c[165]) / (0.9273 - c[165])
    ical_f_hat = np.minimum(1, np.maximum(ical_f_hat_val, 0))

    # -----------------------------------------------------------------
    # IKs
    # -----------------------------------------------------------------
    iksp = y[40]
    # Fraction of phosphorylated IKs
    fp_iks = (iksp + c[144]) / c[143]
    fp_iks = max(0.0001, fp_iks)
    fp_iks = min(0.9999, fp_iks)

    # Effective fraction of phosphorylated IKs channels
    iks_f_hat_val = (fp_iks - c[166]) / (0.785 - c[166])
    iks_f_hat = np.minimum(
        1, np.maximum(iks_f_hat_val, 0)
    )  # np.clip(iks_f_hat_val, 0.0, 1.0)

    # -----------------------------------------------------------------
    # Iup (PLB)
    # -----------------------------------------------------------------
    iup_f_plb = y[41]
    iup_f_pka_val = (iup_f_plb - 0.6662) / (0.9945 - 0.6662)
    iup_f_pka = np.minimum(
        1, np.maximum(iup_f_pka_val, 0)
    )  # np.clip(iup_f_pka_val, 0.0, 1.0)

    # -----------------------------------------------------------------
    # Tni
    # -----------------------------------------------------------------
    f_tni = y[42]
    # Effective fraction of phosphorylated Troponin
    calcium_fhat_val = (f_tni - 0.67352) / (0.99918 - 0.67352)
    calcium_fhat = np.minimum(
        1, np.maximum(calcium_fhat_val, 0)
    )  # np.clip(calcium_fhat_val, 0.0, 1.0)

    # -----------------------------------------------------------------
    # INa
    # -----------------------------------------------------------------
    ina_f_ina = y[43]
    # Effective fraction of phosphorylated INa channels
    ina_f_pka_val = (ina_f_ina - 0.23948) / (0.95014 - 0.23948)
    ina_f_pka = np.minimum(
        1, np.maximum(ina_f_pka_val, 0)
    )  # np.clip(ina_f_pka_val, 0.0, 1.0)

    # -----------------------------------------------------------------
    # INaK
    # -----------------------------------------------------------------
    f_inak = y[44]
    # Effective fraction of phosphorylated INaK pumps
    inak_fhat_val = (f_inak - 0.12635) / (0.99801 - 0.12635)
    inak_fhat = np.minimum(
        1, np.maximum(inak_fhat_val, 0)
    )  # np.clip(inak_fhat_val, 0.0, 1.0)

    # -----------------------------------------------------------------
    # RyR
    # -----------------------------------------------------------------
    ryrp = y[45]
    # Fraction of phosphorylated RyR channels
    fp_ryr = (ryrp + c[160]) / c[150]
    fp_ryr = max(0.0001, fp_ryr)
    fp_ryr = min(0.9999, fp_ryr)

    # Effective fraction of phosphorylated RyR channels
    irel_fhat_val = (fp_ryr - c[164]) / (0.9586 - c[164])
    irel_fhat = np.minimum(
        1, np.maximum(irel_fhat_val, 0)
    )  # np.clip(irel_fhat_val, 0.0, 1.0)

    # -----------------------------------------------------------------
    # IKur
    # -----------------------------------------------------------------
    f_ikur = y[46]
    # Effective fraction of phosphorylated IKur channels
    ikur_fhat_val = (f_ikur - 5.89380e-02) / (0.39375 - 5.89380e-02)
    ikur_fhat = np.minimum(
        1, np.maximum(ikur_fhat_val, 0)
    )  # np.clip(ikur_fhat_val, 0.0, 1.0)

    # -----------------------------------------------------------------
    # Assign final values
    # -----------------------------------------------------------------
    fICaLP = ical_f_hat
    fIKsP = iks_f_hat
    fPLBP = iup_f_pka
    fTnIP = calcium_fhat
    fINaP = ina_f_pka
    fINaKP = inak_fhat
    fRyRP = irel_fhat
    fIKurP = ikur_fhat

    if original_output:
        # Orginal order corresponds to getEffectiveFraction in MATLAB
        output[0] = fICaLP
        output[1] = fIKsP
        output[2] = fPLBP
        output[3] = fTnIP
        output[4] = fINaP
        output[5] = fINaKP
        output[6] = fRyRP  # This is not currently used in the model
        output[7] = fIKurP  # This is not currently used in the model

    else:
        # Order corresponds to PKA_P input in model_TWorld_shock
        output[0] = fINaP
        output[1] = fICaLP
        output[2] = fINaKP
        output[3] = fIKsP
        output[4] = fPLBP
        output[5] = fTnIP
        output[6] = fTnIP

    return output
