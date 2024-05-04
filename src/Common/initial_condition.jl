"""
   Initial condition computation
"""

"""
   Initial profiles and surface values as defined by KiD setup
"""
function init_profile(::Type{FT}, kid_params, thermo_params, z; dry = false) where {FT}

    z_0::FT = kid_params.z_0   # 0.0
    z_1::FT = kid_params.z_1   # 740.0
    z_2::FT = kid_params.z_2   # 3260.0
    rv_0::FT = kid_params.rv_0 # 0.015
    rv_1::FT = kid_params.rv_1 # 0.0138
    rv_2::FT = kid_params.rv_2 # 0.0024
    θ_0::FT = kid_params.tht_0   # 297.9
    θ_1::FT = kid_params.tht_1   # 297.9
    θ_2::FT = kid_params.tht_2   # 312.66

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
    p_0 = kid_params.p0
    SDM_θ_dry_0 = SDM_θ_dry(thermo_params, θ_0, qv_0)
    SDM_ρ_dry_0 = SDM_ρ_dry(thermo_params, p_0, qv_0, θ_0)
    SDM_T_0 = SDM_T(thermo_params, SDM_θ_dry_0, SDM_ρ_dry_0)
    SDM_ρ_0 = SDM_ρ_of_ρ_dry(SDM_ρ_dry_0, qv_0)

    return (qv = qv, θ_std = θ_std, ρ_0 = SDM_ρ_0, z_0 = z_0, z_2 = z_2)
end

"""
    Density derivative as a function of height (assuming no cloud condensate)
    Needed to solve the initial value problem to create the density profile.
"""
function dρ_dz!(ρ, ode_settings, z)

    FT = eltype(ρ)
    (; dry, kid_params, thermo_params) = ode_settings

    # initial profiles
    init = init_profile(FT, kid_params, thermo_params, z, dry = dry)
    θ_std::FT = init.θ_std
    q_vap::FT = init.qv

    # constants
    g::FT = TD.Parameters.grav(thermo_params)
    R_d::FT = TD.Parameters.R_d(thermo_params)
    R_v::FT = TD.Parameters.R_v(thermo_params)
    cp_d::FT = TD.Parameters.cp_d(thermo_params)
    cp_v::FT = TD.Parameters.cp_v(thermo_params)

    θ_dry::FT = SDM_θ_dry(thermo_params, θ_std, q_vap)
    ρ_dry::FT = SDM_ρ_dry_of_ρ(ρ, q_vap)
    T::FT = SDM_T(thermo_params, θ_dry, ρ_dry)
    p::FT = SDM_p(thermo_params, ρ_dry, T, q_vap)

    r_vap::FT = q_vap / (1 - q_vap)
    R_m = R_v / (1 / r_vap + 1) + R_d / (1 + r_vap)
    cp_m = cp_v / (1 / r_vap + 1) + cp_d / (1 + r_vap)

    ρ_SDM = p / R_m / T

    return g / T * ρ_SDM * (R_m / cp_m - 1) / R_m
end

"""
    Solve the initial value problem for the density profile
"""
function ρ_ivp(::Type{FT}, kid_params, thermo_params; dry = false) where {FT}

    init_surface = init_profile(FT, kid_params, thermo_params, 0.0, dry = dry)

    ρ_0::FT = init_surface.ρ_0
    z_0::FT = init_surface.z_0
    z_max::FT = init_surface.z_2
    z_span = (z_0, z_max)

    ode_settings = (; dry, kid_params, thermo_params)

    prob = ODE.ODEProblem(dρ_dz!, ρ_0, z_span, ode_settings)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-10, abstol = 1e-10)

    return sol
end

"""
    Populate the remaining profiles based on the KiD initial condition
    and the density profile
"""
function initial_condition_1d(
    ::Type{FT},
    common_params,
    kid_params,
    thermo_params,
    ρ_profile,
    z;
    dry = false,
) where {FT}

    q_vap::FT = init_profile(FT, kid_params, thermo_params, z, dry = dry).qv
    θ_std::FT = init_profile(FT, kid_params, thermo_params, z, dry = dry).θ_std

    ρ::FT = ρ_profile(z)
    ρ_dry::FT = SDM_ρ_dry_of_ρ(ρ, q_vap)
    ρ_SDP = 1.225

    # assuming no cloud condensate in the initial profile
    θ_liq_ice::FT = θ_std # TODO - compute this based on TS
    q_tot::FT = q_vap
    ρq_tot::FT = q_tot * ρ

    θ_dry::FT = SDM_θ_dry(thermo_params, θ_std, q_vap)
    T::FT = SDM_T(thermo_params, θ_dry, ρ_dry)

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
    ρq_vap::FT = q_vap * ρ
    N_liq::FT = FT(0)
    N_rai::FT = FT(0)
    N_aer_0::FT = common_params.prescribed_Nd * ρ_dry / ρ_SDP
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

    S_Na_activation::FT = FT(0)
    S_Nl_activation::FT = FT(0)

    term_vel_rai::FT = FT(0)
    term_vel_sno::FT = FT(0)
    term_vel_N_rai::FT = FT(0)

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
        ρq_vap,
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
        S_Na_activation,
        S_Nl_activation,
        term_vel_rai,
        term_vel_sno,
        term_vel_N_rai,
    )
end

function initial_condition_1d(
    ::Type{FT},
    ::Val{NM},
    common_params,
    kid_params,
    thermo_params,
    ρ_profile,
    z;
    dry = false,
) where {NM, FT}
    moments::NTuple{NM, FT} = ntuple(_ -> FT(0), NM)
    S_moments::NTuple{NM, FT} = ntuple(_ -> FT(0), NM)
    weighted_vt::NTuple{NM, FT} = ntuple(_ -> FT(0), NM)
    S_ρq_vap::FT = FT(0)
    S_activation::NTuple{NM, FT} = ntuple(_ -> FT(0), NM)
    
    if NM % 3 > 0
        error("must use a multiple of 3 moments")
    end

    # set up the distributions # TODO: make this general for types of moments
    # if NM % 3 == 0 # all gamma
    #     ndist = Int(NM / 3)
    #     nexp = 0
    # elseif NM % 3 == 2 # one exp, remaining gamma
    #     ndist = Int((NM - 2) / 3 + 1)
    #     nexp = 1
    # else # two exp, remaining gamma
    #     ndist = Int((NM - 4) / 3 + 2)
    #     nexp = 2
    # end
    # pdists::NTuple{ND, <: CL.ParticleDistributions.PrimitiveParticleDistribution{FT}} = ntuple(ND) do k
    #     if k <= nexp
    #         CL.ParticleDistributions.ExponentialPrimitiveParticleDistribution(FT(0), FT(k))
    #     else
    #         CL.ParticleDistributions.GammaPrimitiveParticleDistribution(FT(0), FT(k), FT(1))
    #     end
    # end
    pdists::NTuple{Int(NM / 3), CL.ParticleDistributions.GammaPrimitiveParticleDistribution{FT}} = ntuple(Int(NM/3)) do k
        CL.ParticleDistributions.GammaPrimitiveParticleDistribution(FT(0), FT(k), FT(1))
    end
    
    return merge(
        initial_condition_1d(FT, common_params, kid_params, thermo_params, ρ_profile, z, dry=dry),
        (; moments=moments, S_moments=S_moments, S_ρq_vap=S_ρq_vap, S_activation=S_activation, pdists=pdists, weighted_vt=weighted_vt)
    )
end

"""
    Populate the remaining profiles based on given initial conditions including total specific water
    content (liquid + rain) and total number concentration
"""
function initial_condition(::Type{FT}, thermo_params, qt::FT, Nd::FT, k::FT, ρ_dry::FT, z) where {FT}

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
    ρq_ice::FT = FT(0)
    ρq_sno::FT = FT(0)

    ρ = ρ_dry + ρq_liq

    q_liq::FT = ρq_liq / ρ
    q_rai::FT = ρq_rai / ρ
    q_tot::FT = q_liq
    q_ice::FT = FT(0)
    q_sno::FT = FT(0)

    num_ratio = SF.gamma_inc(thrshld, k)[2]
    N_liq::FT = num_ratio * Nd
    N_rai::FT = Nd - N_liq
    N_aer_0::FT = FT(0.0)
    N_aer::FT = FT(0.0)

    T::FT = FT(300)
    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
    p = TD.air_pressure(thermo_params, ts)
    θ_liq_ice = TD.liquid_ice_pottemp(thermo_params, ts)
    θ_dry = TD.dry_pottemp(thermo_params, T, ρ_dry)

    S_ql_moisture::FT = FT(0.0)
    S_qi_moisture::FT = FT(0.0)

    S_qt_precip::FT = FT(0.0)
    S_ql_precip::FT = FT(0.0)
    S_qi_precip::FT = FT(0.0)
    S_qr_precip::FT = FT(0.0)
    S_qs_precip::FT = FT(0.0)
    S_Nl_precip::FT = FT(0)
    S_Nr_precip::FT = FT(0)

    S_Na_activation::FT = FT(0)
    S_Nl_activation::FT = FT(0)

    term_vel_rai::FT = FT(0)
    term_vel_sno::FT = FT(0)
    term_vel_N_rai::FT = FT(0)

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
        S_Na_activation,
        S_Nl_activation,
        term_vel_rai,
        term_vel_sno,
        term_vel_N_rai,
    )
end
