import numpy as np
import numba
import getEffectiveFraction
import getPKASignalling

names_signalling = (
    "fINa_PKA_in",
    "fICaL_PKA_in",
    "fINaK_PKA_in",
    "fIKs_PKA_in",
    "fPLB_PKA_in",
    "fTnI_PKA_in",
    "fMyBPC_PKA_in",
    "Whole_cell_PP1_in",
)

get_pka_signalling = numba.njit(getPKASignalling.get_pka_signalling)
get_effective_fraction = numba.njit(getEffectiveFraction.get_effective_fraction)


@numba.njit
def update_fraction_parameters(
    runSignalingPathway: bool,
    dt: float,
    X0: np.ndarray,
    const_signaling: np.ndarray,
):
    fraction = np.zeros(len(names_signalling))
    if runSignalingPathway:
        # Use forward euler to solve the signaling pathway
        dX_Signaling = get_pka_signalling(X0, const_signaling)
        X0 += dX_Signaling * dt

        get_effective_fraction(y=X0, c=const_signaling, output=fraction[:-1])
        # Concentration of uninhibited PP1 in the cytosolic compartment
        pp1_PP1f_cyt_sum = const_signaling[37] - const_signaling[36] + X0[38]
        # Concentration of uninhibited PP1 in the cytosolic compartment
        PP1f_cyt = 0.5 * (
            np.sqrt(
                pp1_PP1f_cyt_sum**2.0 + 4.0 * const_signaling[37] * const_signaling[36]
            )
            - pp1_PP1f_cyt_sum
        )
        Whole_cell_PP1 = (
            const_signaling[35] / const_signaling[4]
            + const_signaling[34] / const_signaling[5]
            + PP1f_cyt / const_signaling[6]
        )
        fraction[-1] = Whole_cell_PP1
    else:
        get_effective_fraction(y=X0, c=const_signaling, output=fraction[:-1])
        fraction[-1] = 0.13698
    return fraction
