"""
    effective_fractions!(output, u, p)

Compute the 8 effective phosphorylation fractions from a signaling state vector.
Writes `[fICaL, fIKs, fPLB, fTnI, fINa, fINaK, fRyR, fIKur]` into `output`.

# Arguments
- `output`: output vector (8 elements), modified in place
- `u`: signaling state vector (57 elements)
- `p`: parameter vector (167 elements) from [`compute_parameters`](@ref)
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
    ikur_fhat_val = (f_ikur - 5.893798e-2) / (0.393747 - 5.893798e-2)
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
