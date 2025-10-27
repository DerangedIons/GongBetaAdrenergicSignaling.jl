"""
    effective_fractions!(output::AbstractArray, states::AbstractArray, c::AbstractArray)

Compute effective phosphorylation fractions for 8 cardiac substrate proteins.

This function calculates the effective fractions of phosphorylated ion channels and
regulatory proteins that are relevant for cardiac electrophysiology coupling. The
fractions represent the proportion of each substrate in its phosphorylated (active) state
after PKA-mediated phosphorylation.

# Arguments
- `output::AbstractArray`: 8-element array to store computed fractions (modified in-place)
- `states::AbstractArray`: 57-element state vector from the signaling model
- `c::AbstractArray`: 167-element parameter vector (in reference parameter order c1-c167)

# Output Indices
The `output` array contains the following effective fractions:
1. `fICaLP`: L-type calcium channel (ICaL) phosphorylation
2. `fIKsP`: Slow delayed rectifier potassium channel (IKs) phosphorylation
3. `fPLBP`: Phospholamban (PLB) phosphorylation (affects SERCA pump)
4. `fTnIP`: Troponin I (TnI) phosphorylation (affects myofilament Ca²⁺ sensitivity)
5. `fINaP`: Sodium channel (INa) phosphorylation
6. `fINaKP`: Sodium-potassium pump (INaK) phosphorylation
7. `fRyRP`: Ryanodine receptor (RyR) phosphorylation (affects Ca²⁺ release)
8. `fIKurP`: Ultra-rapid delayed rectifier potassium channel (IKur) phosphorylation

# Returns
`nothing` (the function modifies `output` in-place)

# Notes
- Parameters must be in "reference order" (c1-c167) as computed by `compute_parameters`
  or `ConstantsSignalingMyokit2!` from the reference implementation
- States must be in reference order (indices 40-47 contain substrate phosphorylation states)
- All fractions are clamped to valid ranges [0.0, 1.0]
- Based on Gong et al. (2020) beta-adrenergic signaling model

# Reference
Gong, J.Q.X., Susilo, M.E., Sher, A., Musante, C.J., & Sobie, E.A. (2020).
Quantitative analysis of variability in an integrated model of human ventricular
electrophysiology and β-adrenergic signaling. Journal of Molecular and Cellular Cardiology.
DOI: https://doi.org/10.1016/j.yjmcc.2020.04.009
"""
function effective_fractions!(
    output::AbstractArray, states::AbstractArray, c::AbstractArray
)
    # Calculating effective fraction of phosphorylated substrates
    # ICaL
    ICaLp = states[40]
    fp_ICaL_val = (ICaLp + c[163]) / c[156]  # Fraction of phosphorylated ICaL channels
    fp_ICaL = clamp(fp_ICaL_val, 0.0001, 0.9999)
    ical_f_hat_val = (fp_ICaL - c[166]) / (0.9273 - c[166])
    ical_f_hat = clamp(ical_f_hat_val, 0.0, 1.0)  # Effective fraction of phosphorylated ICaL channels

    # IKs
    IKsp = states[41]
    fp_iks_val = (IKsp + c[145]) / c[144]  # Fraction of phosphorylated IKs
    fp_iks = clamp(fp_iks_val, 0.0001, 0.9999)

    iks_f_hat_val = (fp_iks - c[167]) / (0.785 - c[167])
    iks_f_hat = clamp(iks_f_hat_val, 0.0, 1.0)  # Effective fraction of phosphorylated IKs channels

    # Iup (PLB)
    iup_f_plb = states[42]
    #iup_f_pka_val = (iup_f_plb - 0.6591) / (0.9945 - 0.6591)
    iup_f_pka_val = (iup_f_plb - 0.6662) / (0.9945 - 0.6662)
    iup_f_pka = clamp(iup_f_pka_val, 0.0, 1.0)

    # Tni
    f_tni = states[43]
    calcium_fhat_val = (f_tni - 0.6735188) / (0.9991797 - 0.6735188)
    calcium_fhat = clamp(calcium_fhat_val, 0.0, 1.0) # Effective fraction of phosphorylated Troponin

    # INa
    ina_f_ina = states[44]
    ina_f_pka_val = (ina_f_ina - 0.2394795) / (0.9501431 - 0.2394795)
    ina_f_pka = clamp(ina_f_pka_val, 0.0, 1.0) # Effective fraction of phosphorylated INa channels

    # INaK
    f_inak = states[45]
    inak_fhat_val = (f_inak - 0.1263453) / (0.9980137 - 0.1263453)
    inak_fhat = clamp(inak_fhat_val, 0.0, 1.0) # Effective fraction of phosphorylated INaK pumps

    # RyR
    RyRp = states[46]
    fp_RyR_val = (RyRp + c[161]) / c[151]  # Fraction of phosphorylated RyR channels
    fp_RyR = clamp(fp_RyR_val, 0.0001, 0.9999)
    irel_fhat_val = (fp_RyR - c[165]) / (0.9586 - c[165])
    irel_fhat = clamp(irel_fhat_val, 0.0, 1.0)  # Effective fraction of phosphorylated ryr channels

    # IKur
    f_ikur = states[47]
    ikur_fhat_val = (f_ikur - 5.893798e-02) / (0.393747 - 5.893798e-02)
    ikur_fhat = clamp(ikur_fhat_val, 0.0, 1.0) # Effective fraction of phosphorylated IKur channels

    # Some helper to map the names
    fICaLP = ical_f_hat
    fIKsP = iks_f_hat
    fPLBP = iup_f_pka
    fTnIP = calcium_fhat
    fINaP = ina_f_pka
    fINaKP = inak_fhat
    fRyRP = irel_fhat
    fIKurP = ikur_fhat

    output[1] = fICaLP
    output[2] = fIKsP
    output[3] = fPLBP
    output[4] = fTnIP
    output[5] = fINaP
    output[6] = fINaKP
    output[7] = fRyRP
    output[8] = fIKurP

    return nothing
end
