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
    error("initailisation not implemented for a given $(typeof(sm).name.name) and $(typeof(sp).name.name)")
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
        ρq_sno = initial_profiles.ρq_sno,
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
        ρq_sno = initial_profiles.ρq_sno,
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
        q_vap = initial_profiles.q_vap,
        ρq_vap = initial_profiles.ρq_vap,
        q_rai = initial_profiles.q_rai,
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
    # helper function to copy zero field to all fields of a `@NamedTuple{NT}`-valued Field
    zero_nt(zero, ::Type{NT}) where {NT} = NT(ntuple(_ -> copy(zero), fieldcount(NT)))
    zero_nt(::Type{NT}) where {NT} = @. zero_nt(ip.zero, NT)

    # Create a thermo state for aux
    # Allocate the cloud_sources which is a field of containers (tuples)
    # for storing non-equilibrium cloud formation sources
    if moisture isa EquilibriumMoisture
        ts = @. TD.PhaseEquil_ρθq(thermo_params, ip.ρ, ip.θ_liq_ice, ip.q_tot)
        cloud_sources = nothing
    elseif moisture isa NonEquilibriumMoisture || moisture isa MoistureP3
        q = @. TD.PhasePartition(ip.q_tot, ip.q_liq, ip.q_ice)
        ts = @. TD.PhaseNonEquil_ρθq(thermo_params, ip.ρ, ip.θ_liq_ice, q)

        cloud_sources = zero_nt(NamedTuple{(:q_liq, :q_ice)})
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
        tmp = similar(ip.q_tot), tmp2 = similar(ip.q_tot), tmp3 = similar(ip.q_tot),
        tmp_surface = similar(CC.Fields.level(ip.q_tot, 1), Tuple{FT, FT}),
    )

    # Allocate a tuple of fields for storing precomputed thermodynamic variables
    thermo_variables = (; ts, ip.ρ, ip.ρ_dry, ip.p, ip.T, ip.θ_liq_ice, ip.θ_dry)

    microph_variables = (; ip.q_tot, ip.q_liq, ip.q_ice)
    if precip isa Union{NoPrecipitation, Precipitation0M}
        # no additional `microph_variables`
        velocities = nothing
        precip_sources = zero_nt(NamedTuple{(:q_tot, :q_liq, :q_ice)})
        activation_sources = nothing
    elseif precip isa Precipitation1M
        microph_variables = merge(microph_variables, (; ip.q_rai, ip.q_sno))
        velocities = zero_nt(NamedTuple{(:term_vel_rai, :term_vel_sno)})
        precip_sources = zero_nt(NamedTuple{(:q_tot, :q_liq, :q_ice, :q_rai, :q_sno)})
        activation_sources = nothing
    elseif precip isa Precipitation2M
        microph_variables = merge(microph_variables, (; ip.q_rai, ip.q_sno, ip.N_liq, ip.N_rai, ip.N_aer))
        velocities = zero_nt(NamedTuple{(:term_vel_N_rai, :term_vel_rai)})
        precip_sources = zero_nt(NamedTuple{(:q_tot, :q_liq, :q_rai, :N_aer, :N_liq, :N_rai)})
        activation_sources = zero_nt(NamedTuple{(:N_aer, :N_liq)})
    elseif precip isa PrecipitationP3
        microph_variables = (; 
            ip.q_tot, ip.q_liq, ip.q_rai, ip.q_ice, ip.q_rim, ip.q_liqonice, 
            ip.ρq_tot, ip.ρq_liq, ip.ρq_rai, ip.ρq_ice, ip.ρq_rim, ip.ρq_liqonice, 
            ip.B_rim, ip.N_liq, ip.N_rai, ip.N_ice, ip.N_aer, 
            ip.q_vap, ip.ρq_vap,
        )
        velocities = zero_nt(NamedTuple{(:term_vel_rai, :term_vel_ice, :term_vel_N_rai, :term_vel_N_ice)})
        precip_sources = zero_nt(
            NamedTuple{(
                :q_tot, :q_liq, :q_rai, :q_ice, :q_rim, :q_liqonice, :N_aer, :N_liq, :N_rai, :N_ice, :B_rim,
            )},
        )
        p3_boundary_condition_eltype = @NamedTuple{
            ice_start::Bool, _magnitude::FT, _q_flux::FT, _N_flux::FT, _F_rim::FT, _F_liq::FT, _ρ_r_init::FT,
        }
        p3_boundary_condition = @. p3_boundary_condition_eltype(
            tuple(
                copy(precip.p3_boundary_condition.ice_start),
                copy(precip.p3_boundary_condition._magnitude),
                copy(precip.p3_boundary_condition._q_flux),
                copy(precip.p3_boundary_condition._N_flux),
                copy(precip.p3_boundary_condition._F_rim),
                copy(precip.p3_boundary_condition._F_liq),
                copy(precip.p3_boundary_condition._ρ_r_init),
            ),
        )

        activation_sources = nothing
    elseif precip isa CloudyPrecip
        microph_variables = (; 
            ip.q_tot, ip.q_liq, ip.q_ice, ip.q_rai, ip.N_liq, ip.N_rai, ip.N_aer, ip.pdists, ip.moments,
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
    elseif precip isa PrecipitationP3
        aux = merge(aux, (; p3_boundary_condition))
    end

    return aux
end
