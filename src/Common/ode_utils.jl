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
    return NamedTuple{(:ρq_tot,)}.(initial_profiles)
end
function initialise_state(::NonEquilibriumMoisture, ::Union{NoPrecipitation, Precipitation0M}, initial_profiles)
    return NamedTuple{(:ρq_tot, :ρq_liq, :ρq_ice)}.(initial_profiles)
end
function initialise_state(::EquilibriumMoisture, ::Precipitation1M, initial_profiles)
    return NamedTuple{(:ρq_tot, :ρq_rai, :ρq_sno)}.(initial_profiles)
end
function initialise_state(::NonEquilibriumMoisture, ::Precipitation1M, initial_profiles)
    return NamedTuple{(:ρq_tot, :ρq_liq, :ρq_ice, :ρq_rai, :ρq_sno)}.(initial_profiles)
end
function initialise_state(::EquilibriumMoisture, ::Precipitation2M, initial_profiles)
    return NamedTuple{(:ρq_tot, :ρq_rai, :ρq_sno, :N_liq, :N_rai, :N_aer)}.(initial_profiles)
end
function initialise_state(::NonEquilibriumMoisture, ::Precipitation2M, initial_profiles)
    return NamedTuple{(:ρq_tot, :ρq_liq, :ρq_ice, :ρq_rai, :ρq_sno, :N_liq, :N_rai, :N_aer)}.(initial_profiles)
end
function initialise_state(::MoistureP3, ::PrecipitationP3, initial_profiles)
    return NamedTuple{(
        :ρq_tot, :ρq_liq, :ρq_rai, :ρq_ice, :ρq_rim, :ρq_liqonice, :B_rim,
        :N_liq, :N_rai, :N_ice, :N_aer, :q_vap, :ρq_vap, :q_rai,
    )}.(initial_profiles)
end
function initialise_state(::NonEquilibriumMoisture, ::Precipitation2M_P3, initial_profiles)
    warm_state = NamedTuple{(:ρq_tot, :ρq_liq, :ρq_rai, :N_liq, :N_rai, :N_aer)}.(initial_profiles)
    cold_state = NamedTuple{(:ρq_ice, :N_ice, :ρq_rim, :B_rim)}.(initial_profiles)
    tmp = NamedTuple{(:ρq_sno,)}.(initial_profiles)  # only needed to use existing `Precipitation2M` code, should always be zero.
    return merge.(warm_state, cold_state, tmp)
end
function initialise_state(::CloudyMoisture, ::CloudyPrecip, initial_profiles)
    return NamedTuple{(:ρq_vap, :N_aer, :moments)}.(initial_profiles)
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
    zero_nt(zero, ::Type{NT}) where {NT} = NT(ntuple(Returns(zero), Val(fieldcount(NT))))
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
    thermo_ip = NamedTuple{(:ρ, :ρ_dry, :p, :T, :θ_liq_ice, :θ_dry)}.(ip)
    thermo_variables = merge.(NamedTuple{(:ts,)}.(tuple.(ts)), thermo_ip)

    microph_variables = NamedTuple{(:q_tot, :q_liq, :q_ice)}.(ip)
    if precip isa Union{NoPrecipitation, Precipitation0M}
        # no additional `microph_variables`
        velocities = nothing
        precip_sources = zero_nt(NamedTuple{(:q_tot, :q_liq, :q_ice)})
        activation_sources = nothing
    elseif precip isa Precipitation1M
        microph_variables = merge.(microph_variables, NamedTuple{(:q_rai, :q_sno)}.(ip))
        velocities = zero_nt(NamedTuple{(:term_vel_rai, :term_vel_sno)})
        precip_sources = zero_nt(NamedTuple{(:q_tot, :q_liq, :q_ice, :q_rai, :q_sno)})
        activation_sources = nothing
    elseif precip isa Precipitation2M || precip isa Precipitation2M_P3
        microph_variables = merge.(microph_variables, NamedTuple{(:q_rai, :q_sno, :N_liq, :N_rai, :N_aer)}.(ip))
        velocities = zero_nt(NamedTuple{(:term_vel_N_rai, :term_vel_rai)})
        precip_sources = zero_nt(NamedTuple{(:q_tot, :q_liq, :q_rai, :N_aer, :N_liq, :N_rai)})
        activation_sources = zero_nt(NamedTuple{(:N_aer, :N_liq)})
        if precip isa Precipitation2M_P3
            # NOTE: `q_sno` from above is not used in the P3 code, but is needed to use existing `Precipitation2M` code.
            ice_microph_variables = NamedTuple{(:F_rim, :ρ_rim, :logλ)}.(ip)  # variables needed to be precomputed 
            ice_output_variables = NamedTuple{(:N_ice, :ρq_rim, :B_rim)}.(ip)  # variables for netCDF output
            microph_variables = merge.(microph_variables, ice_microph_variables, ice_output_variables)
            ice_velocities = zero_nt(NamedTuple{(:term_vel_N_ice, :term_vel_q_ice)})
            velocities = merge.(velocities, ice_velocities)  # precompute ice terminal velocities
            ice_precip_sources = zero_nt(NamedTuple{(:q_ice, :q_rim, :N_ice, :B_rim)})
            precip_sources = merge.(precip_sources, ice_precip_sources)
        end
    elseif precip isa PrecipitationP3
        microph_variables =
            NamedTuple{(
                :q_tot, :q_liq, :q_rai, :q_ice, :q_rim, :q_liqonice,
                :ρq_tot, :ρq_liq, :ρq_rai, :ρq_ice, :ρq_rim, :ρq_liqonice,
                :B_rim, :N_liq, :N_rai, :N_ice, :N_aer,
                :q_vap, :ρq_vap,
            )}.(ip)
        velocities = zero_nt(NamedTuple{(:term_vel_rai, :term_vel_ice, :term_vel_N_rai, :term_vel_N_ice)})
        precip_sources = zero_nt(
            NamedTuple{(
                :q_tot, :q_liq, :q_rai, :q_ice, :q_rim, :q_liqonice, :N_aer, :N_liq, :N_rai, :N_ice, :B_rim,
            )},
        )
        p3_boundary_condition = precip.p3_boundary_condition
        activation_sources = nothing
    elseif precip isa CloudyPrecip
        microph_variables = NamedTuple{(:q_tot, :q_liq, :q_ice, :q_rai, :N_liq, :N_rai, :N_aer, :pdists, :moments)}.(ip)
        velocities = zero_nt(ip.cloudy_moments_zero, NamedTuple{(:weighted_vt,)})
        precip_sources = NamedTuple{(:moments, :ρq_vap)}.(tuple.(ip.cloudy_moments_zero, ip.zero))
        activation_sources =
            NamedTuple{(:activation, :N_aer, :ρq_vap)}.(tuple.(ip.cloudy_moments_zero, ip.zero, ip.zero))
        cloudy_variables = NamedTuple{(:nm_cloud,)}.(tuple.(Val(cloudy_params.NProgMoms[1])))
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
