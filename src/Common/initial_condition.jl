"""
   Initial condition computation
"""

function read_Jouan_sounding()

    # TODO - fix me later
    long_path = "/Users/climanna/clones/KinematicDriver.jl/src/Common/"
    data = DF.readdlm(long_path*"Jouan_initial_condition.txt", skipstart=1)

    input_T = reverse(data[:, 1])
    #Td = reverse(data[:, 2])
    input_qv = reverse(data[:, 3])
    #Qsat = reverse(data[:, 4])
    input_p = reverse(data[:, 5])
    input_z = reverse(data[:, 10])

    T = IT.linear_interpolation(input_z, input_T; extrapolation_bc = IT.Line())
    qv = IT.linear_interpolation(input_z, input_qv; extrapolation_bc = IT.Line())
    p = IT.linear_interpolation(input_z, input_p; extrapolation_bc = IT.Line())

    return (T = T, qv = qv, p = p)
end

"""
   Initial profiles and surface values as defined by Joauan et al
   TODO - add doi to the paper
"""
function init_profile(::Type{FT}, thermo_params, z) where {FT}

    # TODO - double check
    z_0 = 0
    z_2 = 1.5e3

    sounding = read_Jouan_sounding()
    qv = sounding.qv(z)

    # TODO - dry vs moist pressure
    # TODO - should we also provide the phase partition based on the sounding
    θ_std = TD.dry_pottemp_given_pressure(thermo_params, sounding.T(z), sounding.p(z))

    # density at the surface
    p_0 = sounding.p(z_0)
    qv_0 = sounding.qv(z_0)
    θ_0 = TD.dry_pottemp_given_pressure(thermo_params, sounding.T(z_0), sounding.p(z_0))

    SDM_θ_dry_0 = SDM_θ_dry(thermo_params, θ_0, qv_0)
    SDM_ρ_dry_0 = SDM_ρ_dry(thermo_params, p_0, qv_0, θ_0)
    SDM_T_0 = SDM_T(thermo_params, SDM_θ_dry_0, SDM_ρ_dry_0)
    SDM_ρ_0 = SDM_ρ_of_ρ_dry(SDM_ρ_dry_0, qv_0)

    # TODO - is it ok to assume z0=0?
    return (qv = qv, θ_std = θ_std, ρ_0 = SDM_ρ_0, z_0 = 0, z_2 = z_2)
end

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
    # TODO - make it a user option
    #init = init_profile(FT, kid_params, thermo_params, z, dry = dry)
    init = init_profile(FT, thermo_params, z)

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

    θ_dry::FT = SDM_θ_dry(thermo_params, θ_std, q_vap)
    T::FT = SDM_T(thermo_params, θ_dry, ρ_dry)

    ts = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)
    p::FT = TD.air_pressure(thermo_params, ts)

    N_aer::FT = common_params.prescribed_Nd * ρ_dry / ρ_SDP

    q_liq::FT = FT(0) #TD.liquid_specific_humidity(thermo_params, ts)
    q_ice::FT = FT(0) #TD.ice_specific_humidity(thermo_params, ts)
    q_rai::FT = FT(0)
    q_sno::FT = FT(0)
    q_rim::FT = FT(0)
    B_rim::FT = FT(0)

    ρq_tot::FT = q_tot * ρ
    ρq_liq::FT = q_liq * ρ
    ρq_ice::FT = q_ice * ρ
    ρq_rai::FT = q_rai * ρ
    ρq_sno::FT = q_sno * ρ
    ρq_rim::FT = q_rim * ρ
    ρq_vap::FT = q_vap * ρ

    N_liq::FT = FT(0)
    N_ice::FT = FT(0)
    N_rai::FT = FT(0)

    zero::FT = FT(0)

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
        q_rim,
        B_rim,
        N_liq,
        N_ice,
        N_rai,
        N_aer,
        zero,
    )
end

"""
    Populate the remaining profiles based on given initial conditions including total specific water
    content (liquid + rain) and total number concentration
"""
function initial_condition_0d(
    ::Type{FT},
    thermo_params::TD.Parameters.ThermodynamicsParameters{FT},
    qt::FT,
    Nd::FT,
    k::FT,
    ρ_dry::FT,
) where {FT}

    # qt represents specific water content in cloud and rain. The initialization in PySDM is
    # based on initial gamma distributions. This can lead to initial existence of rain; so here
    # we compute variables for any general initial gamma distributions, by assuming absolute
    # fixed radius threshold of 40 um between cloud droplets and raindrops.
    # Lt : total liquid plus rain content (no vapor; or assuming vapor is contained in dry air density)
    Lt::FT = ρ_dry * qt / (1 - qt)

    rhow = FT(1000)
    radius_th::FT = 40 * 1e-6
    thrshld::FT = 4 / 3 * pi * (radius_th)^3 * rhow / (Lt / Nd / k)

    mass_ratio = SF.gamma_inc(thrshld, k + 1)[2]
    ρq_liq::FT = mass_ratio * Lt
    ρq_rai::FT = Lt - ρq_liq
    ρq_tot::FT = Lt
    ρq_ice::FT = FT(0)
    ρq_sno::FT = FT(0)
    ρq_vap::FT = FT(0)

    ρ = ρ_dry + Lt

    q_liq::FT = ρq_liq / ρ
    q_rai::FT = ρq_rai / ρ
    q_tot::FT = qt
    q_ice::FT = FT(0)
    q_sno::FT = FT(0)

    num_ratio = SF.gamma_inc(thrshld, k)[2]
    N_liq::FT = num_ratio * Nd
    N_rai::FT = Nd - N_liq
    N_aer::FT = FT(0)

    T::FT = FT(300)
    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
    p = TD.air_pressure(thermo_params, ts)
    θ_liq_ice = TD.liquid_ice_pottemp(thermo_params, ts)
    θ_dry = TD.dry_pottemp(thermo_params, T, ρ_dry)

    zero::FT = FT(0)

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
        zero,
    )
end

function cloudy_initial_condition(pdists, ip, k = 1)

    NM = sum(CL.ParticleDistributions.nparams.(pdists))

    FT = eltype(ip.ρq_liq)
    Lt::FT = ip.ρq_liq + ip.ρq_rai
    Nd::FT = ip.N_liq + ip.N_rai

    moments::NTuple{NM, FT} = ntuple(NM) do j
        if j == 1
            Nd
        elseif j == 2
            Lt
        elseif j == 3 && CL.ParticleDistributions.nparams(pdists[1]) == 3
            ifelse(Nd < eps(FT), FT(0), Lt^2 / Nd * (k + 1) / k)
        else
            FT(0)
        end
    end

    cloudy_moments_zero::NTuple{NM, FT} = ntuple(_ -> FT(0), NM)

    return merge(ip, (; moments = moments, pdists = pdists, cloudy_moments_zero))
end

function p3_initial_condition(
    ::Type{FT},
    kid_params,
    thermo_params,
    z;
    _q_init = FT(0),
    _N_init = FT(0),
    _F_rim = FT(0),
    _F_liq = FT(0),
    _ρ_r = FT(0),
    z_top,
    ice_start,
    dry = false,
) where {FT}

    # initialize water vapor profile
    # using same setup as initial_condition_1d
    q_vap::FT = init_profile(FT, kid_params, thermo_params, z, dry = dry).qv

    # if ice_start, start with small ice bump
    # otherwise, introduce particles solely through
    # the boundary in advection_tendency!()
    _z_band = z_top * FT(0.2)
    has_ice(z) = ifelse(ice_start, z > (z_top - _z_band), false)
    ice_bump(χ, z) = χ * cos(z - (z_top - 0.5 * _z_band))

    # P3 state variables (minus B_rim -- computed below)
    # (initialized as kg/kg before converting to kg/m3)
    q_ice::FT = has_ice(z) ? ice_bump(_q_init, z) : FT(0)
    N_ice_kg::FT = has_ice(z) ? ice_bump(_N_init, z) : FT(0)
    q_liqonice::FT = _F_liq * q_ice
    q_rim::FT = _F_rim * (q_ice - q_liqonice)

    # 2M warm rain state variables
    q_liq::FT = FT(0)
    N_liq::FT = FT(0)
    q_rai::FT = FT(0)
    N_rai::FT = FT(0)

    # q_tot - single p3 category + 2M categories + vapor
    q_tot::FT = q_ice + q_liq + q_rai + q_vap

    # thermodynamics:
    # for Case 1 of Cholette et al 2019 paper
    # use approximate temperature values from
    # sounding: kin1d/soundings/Tsfc2.txt
    # (path in p3 fortran github repo)
    T::FT = -0.004 * (z - 500) + 273.15 # temperature
    p::FT = 990 - (0.1 * z) # pressure
    _q = TD.PhasePartition(q_tot, q_liq, q_ice)
    ρ::FT = TD.air_density(thermo_params, T, p, _q) # total density
    ρ_dry::FT = TD.air_density(thermo_params, T, p)
    ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, _q) # thermodynamic state

    # P3 scheme uses L (kg/m3) = ρ * q
    # so we compute ρq and pass that to P3
    # also we convert N to (1/m3)
    # since paper init values are given
    # in (kg/kg), (1/kg)
    ρq_tot::FT = ρ * q_tot
    N_ice::FT = ρ * N_ice_kg
    ρq_ice::FT = ρ * q_ice
    ρq_rim::FT = ρ * q_rim
    ρq_liqonice::FT = ρ * q_liqonice
    ρq_rai::FT = ρ * q_rai
    ρq_liq::FT = ρ * q_liq
    ρq_vap::FT = ρ * q_vap

    # also compute B_rim from L_rim, ρ_r
    B_rim::FT = ρq_rim / _ρ_r


    # unused quantities:
    q_sno::FT = FT(0)
    ρq_sno::FT = FT(0)
    θ_liq_ice::FT = TD.liquid_ice_pottemp(thermo_params, ts)
    N_aer::FT = FT(0)
    θ_dry::FT = TD.dry_pottemp(thermo_params, T, ρ)

    # zero:
    zero::FT = FT(0)

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
        ρq_rim,
        ρq_sno,
        ρq_vap,
        ρq_liqonice,
        q_tot,
        q_liq,
        q_ice,
        q_rai,
        q_sno,
        q_rim,
        q_liqonice,
        q_vap,
        B_rim,
        N_liq,
        N_ice,
        N_rai,
        N_aer,
        zero,
    )
end
