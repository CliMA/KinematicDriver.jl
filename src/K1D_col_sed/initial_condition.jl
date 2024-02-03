"""
   Initial condition computation
"""

"""
    Populate the remaining profiles based on the given initial condition
"""
function init_1d_column(::Type{FT}, qt::FT, Nd::FT, k::FT, ρ_dry::FT, z) where {FT}

    # qt represents specific water content in cloud and rain. The initialization in PySDM is 
    # based on initial gamma distributions. This can lead to initial existence of rain; so here
    # we compute variables for any general initial gamma distributions, by assuming absolute
    # fixed radius threshold of 40 um between cloud droplets and raindrops.
    # Ltr : total liquid plus rain content (no vapor; or assuming vapor is contained in dry air density)
    L_tr::FT = ρ_dry * qt / (1 - qt)

    rhow = FT(1000)
    radius_th::FT = 40 * 1e-6
    thrshld::FT = 4 / 3 * pi * (radius_th)^3 * rhow / (L_tr / Nd / k)
    mass_ratio = SF.gamma_inc(thrshld, k + 1)[2]

    ρq_liq::FT = mass_ratio * L_tr
    ρq_rai::FT = L_tr - ρq_liq
    ρq_tot::FT = ρq_liq
    ρ = ρ_dry + ρq_liq

    q_liq::FT = ρq_liq / ρ
    q_rai::FT = ρq_rai / ρ
    q_tot::FT = q_liq

    num_ratio = SF.gamma_inc(thrshld, k)[2]
    N_liq::FT = num_ratio * Nd
    N_rai::FT = Nd - N_liq

    S_qt::FT = FT(0.0)
    S_ql::FT = FT(0.0)
    S_qr::FT = FT(0.0)
    S_Nl::FT = FT(0)
    S_Nr::FT = FT(0)

    return (; ρ, ρ_dry, ρq_tot, ρq_liq, ρq_rai, q_tot, q_liq, q_rai, N_liq, N_rai, S_qt, S_ql, S_qr, S_Nl, S_Nr)
end
