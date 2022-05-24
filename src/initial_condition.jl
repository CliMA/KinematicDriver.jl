"""
   Initial condition computation
"""

"""
   Initial profiles and surface values as defined by KiD setup
"""
function init_condition(::Type{FT}, params, z) where {FT}

    z_0::FT = 0.0
    z_1::FT = 740.0
    z_2::FT = 3260.0
    qv_0::FT = 0.015
    qv_1::FT = 0.0138
    qv_2::FT = 0.0024
    θ_0::FT = 297.9
    θ_1::FT = 297.9
    θ_2::FT = 312.66

    # profile of water vapour specific humidity (TODO - or is it mixing ratio?)
    qv::FT = z < z_1 ? qv_0 + (qv_1 - qv_0)/(z_1 - z_0) * (z - z_0) : qv_1 + (qv_2 - qv_1)/(z_2 - z_1) * (z - z_1)

    # profile of potential temperature
    θ::FT = z < z_1 ? θ_0 : θ_1 + (θ_2 - θ_1)/(z_2 - z_1) * (z - z_1)

    # density at the surface
    p_0::FT = 100200.0
    q_0 = TD.PhasePartition(qv_0, 0.0, 0.0)
    T_0::FT = θ_0 * TD.exner_given_pressure(params, p_0,  q_0)
    ρ_0::FT = TD.air_density(params, T_0, p_0, q_0)

    return(qv = qv, θ = θ, ρ_0 = ρ_0, z_0 = z_0, z_2 = z_2)
end

"""
    TODO
"""
function init_condition_dry(::Type{FT}, params, z) where {FT}

    z_0::FT = 0.0
    z_1::FT = 740.0
    z_2::FT = 3260.0
    qv_0::FT = 0.0
    qv_1::FT = 0.0
    qv_2::FT = 0.0
    θ_0::FT = 297.9
    θ_1::FT = 297.9
    θ_2::FT = 312.66

    # profile of water vapour specific humidity (TODO - or is it mixing ratio?)
    qv::FT = z < z_1 ? qv_0 + (qv_1 - qv_0)/(z_1 - z_0) * (z - z_0) : qv_1 + (qv_2 - qv_1)/(z_2 - z_1) * (z - z_1)

    # profile of potential temperature
    θ::FT = z < z_1 ? θ_0 : θ_1 + (θ_2 - θ_1)/(z_2 - z_1) * (z - z_1)

    # density at the surface
    p_0::FT = 100200.0
    q_0 = TD.PhasePartition(qv_0, 0.0, 0.0)
    T_0::FT = θ_0 * TD.exner_given_pressure(params, p_0,  q_0)
    ρ_0::FT = TD.air_density(params, T_0, p_0, q_0)

    return(qv = qv, θ = θ, ρ_0 = ρ_0, z_0 = z_0, z_2 = z_2)
end

"""
    Density derivative as a function of height (assuming no cloud condensate)
    Needed to solve the initial value problem to create the density profile.
"""
function dρ_dz!(ρ, params, z)

    FT = eltype(ρ)

    # initial profiles
    init = init_condition(FT, params, z)
    θ::FT = init.θ
    qv::FT = init.qv

    q = TD.PhasePartition(qv, 0.0, 0.0)

    # constants
    g::FT = CP.Planet.grav(params)
    cp_m::FT = TD.cp_m(params, q)
    R_m::FT = TD.gas_constant_air(params, q)

    T::FT = θ * (ρ * θ / CP.Planet.MSLP(params) * R_m)^((R_m / cp_m) / (1 - R_m / cp_m))

    return g / T * ρ * (R_m / cp_m - 1) / R_m
end

"""
    Solve the initial value problem for the density profile
"""
function ρ_ivp(::Type{FT}, params; dry=false) where {FT}
    if dry
        init_surface = init_condition(FT, params, 0.0)
    else
        init_surface = init_condition_dry(FT, params, 0.0)
    end

    ρ_0::FT = init_surface.ρ_0
    z_0::FT = init_surface.z_0
    z_max::FT = init_surface.z_2

    z_span = (z_0, z_max)
    prob = ODE.ODEProblem(dρ_dz!, ρ_0, z_span, params)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol=1e-8, abstol=1e-8)

    return sol
end

"""
    Populate the remaining profiles based on the KiD initial condition
    and the density profile
"""
function init_1d_column(::Type{FT}, params, ρ_profile, z; dry=false) where {FT}

    if dry
        q_tot::FT = init_condition_dry(FT, params, z).qv
        θ_liq_ice::FT = init_condition_dry(FT, params, z).θ
    else
        q_tot::FT = init_condition(FT, params, z).qv
        θ_liq_ice::FT = init_condition(FT, params, z).θ
    end

    ρ::FT = ρ_profile(z)
    ρq_tot::FT = q_tot * ρ

    ts = TD.PhaseEquil_ρθq(params, ρ, θ_liq_ice, q_tot)

    T::FT = TD.air_temperature(params, ts)
    θ_dry::FT = TD.dry_pottemp(params, ts)
    p::FT = TD.air_pressure(params, ts)

    q_liq::FT = TD.liquid_specific_humidity(params, ts)
    q_ice::FT = TD.ice_specific_humidity(params, ts)
    ρq_liq::FT = q_liq * ρ
    ρq_ice::FT = q_ice * ρ

    q_rai::FT = FT(0.0)
    q_sno::FT = FT(0.0)
    ρq_rai::FT = q_rai * ρ
    ρq_sno::FT = q_sno * ρ

    return(; ρ, T, p, θ_liq_ice, θ_dry, ρq_tot, ρq_liq, ρq_ice, ρq_rai, ρq_sno, q_tot, q_liq, q_ice, q_rai, q_sno)
end
