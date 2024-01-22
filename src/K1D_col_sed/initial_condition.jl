"""
   Initial condition computation
"""

"""
    Populate the remaining profiles based on the given initial condition
"""
function init_1d_column(::Type{FT}, qt::FT, Nd::FT, ρ_dry::FT, z) where {FT}

    ρ = ρ_dry / (1-qt)
    
    q_tot::FT = qt
    q_liq::FT = qt
    q_rai::FT = FT(0)
    ρq_tot::FT = q_tot * ρ
    ρq_liq::FT = q_liq * ρ
    ρq_rai::FT = q_rai * ρ

    N_liq::FT = Nd
    N_rai::FT = FT(0)

    S_qt::FT = FT(0.0)
    S_ql::FT = FT(0.0)
    S_qr::FT = FT(0.0)
    S_Nl::FT = FT(0)
    S_Nr::FT = FT(0)

    return (;
        ρ,
        ρ_dry,
        ρq_tot,
        ρq_liq,
        ρq_rai,
        q_tot,
        q_liq,
        q_rai,
        N_liq,
        N_rai,
        S_qt,
        S_ql,
        S_qr,
        S_Nl,
        S_Nr,
    )
end