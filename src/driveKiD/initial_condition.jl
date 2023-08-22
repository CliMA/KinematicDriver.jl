"""
   Initial condition computation
"""

"""
   Initial profiles and surface values as defined by KiD setup
"""
function init_condition(::Type{FT}, params, z; dry = false) where {FT}

    z_0::FT = 0.0
    z_1::FT = 740.0
    z_2::FT = 3260.0
    rv_0::FT = 0.015
    rv_1::FT = 0.0138
    rv_2::FT = 0.0024
    θ_0::FT = 297.9
    θ_1::FT = 297.9
    θ_2::FT = 312.66

    if dry
        rv_0 = 0.0
        rv_1 = 0.0
        rv_2 = 0.0
    end

    # profile of water vapour mixing ratio
    rv::FT = z < z_1 ? rv_0 + (rv_1 - rv_0) / (z_1 - z_0) * (z - z_0) : rv_1 + (rv_2 - rv_1) / (z_2 - z_1) * (z - z_1)
    qv::FT = rv / (1 + rv)
    qv_0::FT = rv_0 / (1 + rv_0)

    # profile of potential temperature
    θ_std::FT = z < z_1 ? θ_0 : θ_1 + (θ_2 - θ_1) / (z_2 - z_1) * (z - z_1)

    # density at the surface
    p_0 = params.p0
    SDM_θ_dry_0 = SDM_θ_dry(params, θ_0, qv_0)
    SDM_ρ_dry_0 = SDM_ρ_dry(params, p_0, qv_0, θ_0)
    SDM_T_0 = SDM_T(params, SDM_θ_dry_0, SDM_ρ_dry_0)
    SDM_ρ_0 = SDM_ρ_of_ρ_dry(SDM_ρ_dry_0, qv_0)

    return (qv = qv, θ_std = θ_std, ρ_0 = SDM_ρ_0, z_0 = z_0, z_2 = z_2)
end

"""
    Density derivative as a function of height (assuming no cloud condensate)
    Needed to solve the initial value problem to create the density profile.
"""
function dρ_dz!(ρ, ode_settings, z)

    FT = eltype(ρ)

    dry = ode_settings.dry
    params = ode_settings.params

    # initial profiles
    init = init_condition(FT, params, z, dry = dry)
    θ_std::FT = init.θ_std
    q_vap::FT = init.qv

    # constants
    g::FT = KP.grav(params)
    R_d::FT = KP.R_d(params)
    R_v::FT = KP.R_v(params)
    cp_d::FT = KP.cp_d(params)
    cp_v::FT = KP.cp_v(params)

    θ_dry::FT = SDM_θ_dry(params, θ_std, q_vap)
    ρ_dry::FT = SDM_ρ_dry_of_ρ(ρ, q_vap)
    T::FT = SDM_T(params, θ_dry, ρ_dry)
    p::FT = SDM_p(params, ρ_dry, T, q_vap)

    r_vap::FT = q_vap / (1 - q_vap)
    R_m = R_v / (1 / r_vap + 1) + R_d / (1 + r_vap)
    cp_m = cp_v / (1 / r_vap + 1) + cp_d / (1 + r_vap)

    ρ_SDM = p / R_m / T

    return g / T * ρ_SDM * (R_m / cp_m - 1) / R_m
end

"""
    Solve the initial value problem for the density profile
"""
function ρ_ivp(::Type{FT}, params; dry = false) where {FT}

    init_surface = init_condition(FT, params, 0.0, dry = dry)

    ρ_0::FT = init_surface.ρ_0
    z_0::FT = init_surface.z_0
    z_max::FT = init_surface.z_2
    z_span = (z_0, z_max)

    ode_settings = (; dry, params)

    prob = ODE.ODEProblem(dρ_dz!, ρ_0, z_span, ode_settings)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-10, abstol = 1e-10)

    return sol
end

"""
    Populate the remaining profiles based on the KiD initial condition
    and the density profile
"""
function init_1d_column(::Type{FT}, params, ρ_profile, z; dry = false) where {FT}

    thermo_params = KP.thermodynamics_params(params)

    q_vap::FT = init_condition(FT, params, z, dry = dry).qv
    θ_std::FT = init_condition(FT, params, z, dry = dry).θ_std

    ρ::FT = ρ_profile(z)
    ρ_dry::FT = SDM_ρ_dry_of_ρ(ρ, q_vap)
    ρ_SDP = 1.225

    # assuming no cloud condensate in the initial profile
    θ_liq_ice::FT = θ_std # TODO - compute this based on TS
    q_tot::FT = q_vap
    ρq_tot::FT = q_tot * ρ

    θ_dry::FT = SDM_θ_dry(params, θ_std, q_vap)
    T::FT = SDM_T(params, θ_dry, ρ_dry)

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
    N_aer_0::FT = Parameters.prescribed_Nd(params) * ρ_dry / ρ_SDP
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
