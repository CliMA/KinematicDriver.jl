"""
  Main buliding blocks for the ODE solver
"""

"""
   Struct for storing the:
    - timestepping timestep `dt`,
    - output timestep `dt_io` and
    - simulation time `t_max`.
"""

mutable struct TimeStepping{FT <: Real}
    dt::FT
    dt_io::FT
    t_max::FT
end

"""
    Interface to ODE solver. Returns the ClimaCore.jl FieldVector with
    ODE solver state variables. The state is created via dispatching
    on different moisture and precipitation types
"""
function initialise_state(sm::AbstractMoistureStyle, sp::AbstractPrecipitationStyle, initial_profiles)
    error("initailisation not implemented for a given $sm and $sp")
end
function initialise_state(::EquilibriumMoisture, ::Union{NoPrecipitation, Precipitation0M}, initial_profiles)
    return CC.Fields.FieldVector(; ρq_tot = initial_profiles.ρq_tot)
end
function initialise_state(::NonEquilibriumMoisture, ::Union{NoPrecipitation, Precipitation0M}, initial_profiles)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_liq = initial_profiles.ρq_liq,
        ρq_ice = initial_profiles.ρq_ice,
    )
end
function initialise_state(::EquilibriumMoisture, ::Precipitation1M, initial_profiles)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_rai = initial_profiles.ρq_rai,
        ρq_sno = initial_profiles.ρq_sno,
    )
end
function initialise_state(::NonEquilibriumMoisture, ::Precipitation1M, initial_profiles)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_liq = initial_profiles.ρq_liq,
        ρq_ice = initial_profiles.ρq_ice,
        ρq_rai = initial_profiles.ρq_rai,
        ρq_sno = initial_profiles.ρq_sno,
    )
end
function initialise_state(::EquilibriumMoisture, ::Precipitation2M, initial_profiles)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_rai = initial_profiles.ρq_rai,
        N_liq = initial_profiles.N_liq,
        N_rai = initial_profiles.N_rai,
        N_aer = initial_profiles.N_aer,
    )
end
function initialise_state(::NonEquilibriumMoisture, ::Precipitation2M, initial_profiles)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_liq = initial_profiles.ρq_liq,
        ρq_ice = initial_profiles.ρq_ice,
        ρq_rai = initial_profiles.ρq_rai,
        N_liq = initial_profiles.N_liq,
        N_rai = initial_profiles.N_rai,
        N_aer = initial_profiles.N_aer,
    )
end
function initialise_state(::MoistureP3, ::PrecipitationP3, initial_profiles)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_liq = initial_profiles.ρq_liq,
        ρq_rai = initial_profiles.ρq_rai,
        ρq_ice = initial_profiles.ρq_ice,
        ρq_rim = initial_profiles.ρq_rim,
        ρq_liqonice = initial_profiles.ρq_liqonice,
        B_rim = initial_profiles.B_rim,
        N_liq = initial_profiles.N_liq,
        N_rai = initial_profiles.N_rai,
        N_ice = initial_profiles.N_ice,
        N_aer = initial_profiles.N_aer,
    )
end
function initialise_state(::CloudyMoisture, ::CloudyPrecip, initial_profiles)
    return CC.Fields.FieldVector(;
        ρq_vap = initial_profiles.ρq_vap,
        N_aer = initial_profiles.N_aer,
        moments = initial_profiles.moments,
    )
end

"""
   Interface to ODE solver. It initializes the auxiliary state.
   The auxiliary state is created as a ClimaCore FieldVector
   and passed to ODE solver via the `p` parameter of the ODEProblem.
"""
function initialise_aux(
    FT,
    ip,
    common_params,
    thermo_params,
    air_params,
    activation_params,
    TS,
    Stats,
    moisture,
    precip,
    cloudy_params = nothing,
)

    # Create a thermo state for aux
    # Allocate the cloud_sources which is a field of containers (tuples)
    # for storing non-equilibrium cloud formation sources
    if moisture isa EquilibriumMoisture
        ts = @. TD.PhaseEquil_ρθq(thermo_params, ip.ρ, ip.θ_liq_ice, ip.q_tot)
        cloud_sources = nothing
    elseif moisture isa NonEquilibriumMoisture || moisture isa MoistureP3
        q = @. TD.PhasePartition(ip.q_tot, ip.q_liq, ip.q_ice)
        ts = @. TD.PhaseNonEquil_ρθq(thermo_params, ip.ρ, ip.θ_liq_ice, q)

        cloud_sources_eltype = @NamedTuple{q_liq::FT, q_ice::FT}
        cloud_sources = @. cloud_sources_eltype(tuple(copy(ip.zero), copy(ip.zero)))
    elseif moisture isa CloudyMoisture
        q = @. TD.PhasePartition(ip.q_tot, ip.q_liq, ip.q_ice)
        ts = @. TD.PhaseNonEquil_ρTq(thermo_params, ip.ρ, ip.T, q)
        cloud_sources = nothing
    else
        error(
            "Wrong moisture choice $moisture. The supported options are EquilibriumMoisture, NonEquilibriumMoisture, and CloudyMoisture",
        )
    end

    # Allocate scratch which is a tuple of fields for storing intermediate outputs
    scratch = (;
        tmp = similar(ip.q_tot),
        tmp2 = similar(ip.q_tot),
        tmp3 = similar(ip.q_tot),
        tmp_surface = similar(CC.Fields.level(ip.q_tot, 1), Tuple{FT, FT}),
    )

    # Allocate a tuple of fields for storing precomputed thermodynamic variables
    thermo_variables =
        (; ts = ts, ρ = ip.ρ, ρ_dry = ip.ρ_dry, p = ip.p, T = ip.T, θ_liq_ice = ip.θ_liq_ice, θ_dry = ip.θ_dry)

    if precip isa Union{NoPrecipitation, Precipitation0M}
        microph_variables = (; q_tot = ip.q_tot, q_liq = ip.q_liq, q_ice = ip.q_ice)
        velocities = nothing
        precip_sources_eltype = @NamedTuple{q_tot::FT, q_liq::FT, q_ice::FT}
        precip_sources = @. precip_sources_eltype(tuple(copy(ip.zero), copy(ip.zero), copy(ip.zero)))
        activation_sources = nothing
    elseif precip isa Precipitation1M
        microph_variables = (; q_tot = ip.q_tot, q_liq = ip.q_liq, q_ice = ip.q_ice, q_rai = ip.q_rai, q_sno = ip.q_sno)
        velocities = (; term_vel_rai = copy(ip.zero), term_vel_sno = copy(ip.zero))
        precip_sources_eltype = @NamedTuple{q_tot::FT, q_liq::FT, q_ice::FT, q_rai::FT, q_sno::FT}
        precip_sources =
            @. precip_sources_eltype(tuple(copy(ip.zero), copy(ip.zero), copy(ip.zero), copy(ip.zero), copy(ip.zero)))
        activation_sources = nothing
    elseif precip isa Precipitation2M
        microph_variables = (;
            q_tot = ip.q_tot,
            q_liq = ip.q_liq,
            q_ice = ip.q_ice,
            q_rai = ip.q_rai,
            N_liq = ip.N_liq,
            N_rai = ip.N_rai,
            N_aer = ip.N_aer,
        )
        velocities = (; term_vel_N_rai = copy(ip.zero), term_vel_rai = copy(ip.zero))
        precip_sources_eltype = @NamedTuple{q_tot::FT, q_liq::FT, q_rai::FT, N_aer::FT, N_liq::FT, N_rai::FT}
        precip_sources = @. precip_sources_eltype(
            tuple(copy(ip.zero), copy(ip.zero), copy(ip.zero), copy(ip.zero), copy(ip.zero), copy(ip.zero)),
        )
        activation_sources_eltype = @NamedTuple{N_aer::FT, N_liq::FT}
        activation_sources = @. activation_sources_eltype(tuple(copy(ip.zero), copy(ip.zero)))
    elseif precip isa PrecipitationP3
        microph_variables = (;
            q_tot = ip.q_tot,
            q_liq = ip.q_liq,
            q_rai = ip.q_rai,
            q_ice = ip.q_ice,
            q_rim = ip.q_rim,
            q_liqonice = ip.q_liqonice,
            B_rim = ip.B_rim,
            N_liq = ip.N_liq,
            N_rai = ip.N_rai,
            N_ice = ip.N_ice,
            N_aer = ip.N_aer,
            ρq_tot = ip.ρq_tot,
            ρq_liq = ip.ρq_liq,
            ρq_rai = ip.ρq_rai,
            ρq_ice = ip.ρq_ice,
            ρq_rim = ip.ρq_rim,
            ρq_liqonice = ip.ρq_liqonice,
        )
        velocities = (;
            term_vel_rai = copy(ip.zero),
            term_vel_ice = copy(ip.zero),
            term_vel_N_rai = copy(ip.zero),
            term_vel_N_ice = copy(ip.zero),
        )
        precip_sources_eltype = @NamedTuple{
            q_tot::FT,
            q_liq::FT,
            q_rai::FT,
            q_ice::FT,
            q_rim::FT,
            q_liqonice::FT,
            N_aer::FT,
            N_liq::FT,
            N_rai::FT,
            N_ice::FT,
            B_rim::FT,
        }
        precip_sources = @. precip_sources_eltype(
            tuple(
                copy(ip.zero),
                copy(ip.zero),
                copy(ip.zero),
                copy(ip.zero),
                copy(ip.zero),
                copy(ip.zero),
                copy(ip.zero),
                copy(ip.zero),
                copy(ip.zero),
                copy(ip.zero),
                copy(ip.zero),
            ),
        )

        activation_sources = nothing
    elseif precip isa CloudyPrecip
        microph_variables = (;
            q_tot = ip.q_tot,
            q_liq = ip.q_liq,
            q_ice = ip.q_ice,
            q_rai = ip.q_rai,
            N_liq = ip.N_liq,
            N_rai = ip.N_rai,
            N_aer = ip.N_aer,
            pdists = ip.pdists,
            moments = ip.moments,
        )
        velocities = (; weighted_vt = copy(ip.cloudy_moments_zero))
        precip_sources = (; moments = copy(ip.cloudy_moments_zero), ρq_vap = copy(ip.zero))
        activation_sources =
            (; activation = copy(ip.cloudy_moments_zero), N_aer = copy(ip.zero), ρq_vap = copy(ip.zero))
        cloudy_variables = (; nm_cloud = Val(cloudy_params.NProgMoms[1]))
        scratch = merge(scratch, (; tmp_cloudy = similar(ip.cloudy_moments_zero)))
    else
        error("Wrong precipitation choice $precip")
    end

    aux = (;
        scratch,
        thermo_variables,
        microph_variables,
        activation_sources,
        cloud_sources,
        precip_sources,
        velocities,
        common_params,
        thermo_params,
        air_params,
        activation_params,
        Stats,
        TS,
    )

    if precip isa CloudyPrecip
        aux = merge(aux, (; cloudy_params, cloudy_variables))
    end

    return aux
end
