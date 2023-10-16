"""
   Initial condition computation
"""

"""
    Populate the remaining profiles based on the KiD initial condition
    and the density profile
"""
function init_2d_domain(::Type{FT}, kid_params, thermo_params, ρ_profile, x, z; dry = false) where {FT}

    q_vap::FT = K1D.init_condition(FT, kid_params, thermo_params, z, dry = dry).qv
    θ_std::FT = K1D.init_condition(FT, kid_params, thermo_params, z, dry = dry).θ_std

    ρ::FT = ρ_profile(z)
    ρ_dry::FT = K1D.SDM_ρ_dry_of_ρ(ρ, q_vap)
    ρ_SDP::FT = 1.225

    # assuming no cloud condensate in the initial profile
    θ_liq_ice::FT = θ_std
    q_tot::FT = q_vap
    ρq_tot::FT = q_tot * ρ

    θ_dry::FT = K1D.SDM_θ_dry(thermo_params, θ_std, q_vap)
    T::FT = K1D.SDM_T(thermo_params, θ_dry, ρ_dry)

    ts = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)
    p::FT = TD.air_pressure(thermo_params, ts)

    q_liq::FT = TD.liquid_specific_humidity(thermo_params, ts)
    q_ice::FT = TD.ice_specific_humidity(thermo_params, ts)
    ρq_liq::FT = q_liq * ρ
    ρq_ice::FT = q_ice * ρ

    q_rai::FT = FT(0.0)
    q_sno::FT = FT(0.0)
    ρq_rai::FT = q_rai * ρ
    ρq_sno::FT = q_sno * ρ
    N_liq::FT = FT(0)
    N_rai::FT = FT(0)
    N_aer_0::FT = kid_params.prescribed_Nd * ρ_dry / ρ_SDP
    N_aer::FT = N_aer_0

    S_ql_moisture::FT = FT(0.0)
    S_qi_moisture::FT = FT(0.0)

    S_qt_precip::FT = FT(0.0)
    S_ql_precip::FT = FT(0.0)
    S_qi_precip::FT = FT(0.0)
    S_qr_precip::FT = FT(0.0)
    S_qs_precip::FT = FT(0.0)
    S_Nl_precip::FT = FT(0)
    S_Nr_precip::FT = FT(0)
    S_Na::FT = FT(0)

    return (;
        ρ,
        ρ_dry,
        T,
        p,
        θ_liq_ice,
        θ_dry,
        ρq_tot,
        ρq_liq,
        ρq_ice,
        ρq_rai,
        ρq_sno,
        q_tot,
        q_liq,
        q_ice,
        q_rai,
        q_sno,
        N_liq,
        N_rai,
        N_aer,
        N_aer_0,
        S_ql_moisture,
        S_qi_moisture,
        S_qt_precip,
        S_ql_precip,
        S_qi_precip,
        S_qr_precip,
        S_qs_precip,
        S_Nl_precip,
        S_Nr_precip,
        S_Na,
    )
end
