"""
  Main building blocks of the KiD model
"""

"""
   Interface to ClimaCore.jl Returns an instance of space and face space
   (which are a discretized function space over our computational domain).
"""
function make_function_space(FT, z_min, z_max, n_elem)
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(z_min),
        CC.Geometry.ZPoint{FT}(z_max),
        boundary_tags = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems = n_elem)

    space = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    face_space = CC.Spaces.FaceFiniteDifferenceSpace(space)

    return space, face_space
end

"""
   Interface to ODE solver. Returns the function needed to compute the
   right hand side of the solved ODE. The rhs is assembled via dispatch
   based on the moisture and precipitation types.
"""
function make_rhs_function(ms::AbstractMoistureStyle, ps::AbstractPrecipitationStyle)
    function rhs!(dY, Y, aux, t)

        for eq_style in [ms, ps]
            zero_tendencies!(eq_style, dY, Y, aux, t)
        end

        precompute_aux_prescribed_velocity!(aux, t)
        precompute_aux_thermo!(ms, dY, Y, aux, t)
        precompute_aux_precip!(ps, dY, Y, aux, t)

        for eq_style in [ms, ps]
            advection_tendency!(eq_style, dY, Y, aux, t)
            sources_tendency!(eq_style, dY, Y, aux, t)
        end

    end
    return rhs!
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

"""
   Interface to ODE solver. It initializes the auxiliary state.
   The auxiliary state is created as a ClimaCore FieldVector
   and passed to ODE solver via the `p` parameter of the ODEProblem.
"""
function initialise_aux(FT, ip, params, w_params, TS, Stats, face_space, moisture)
    q_surf = init_condition(FT, params, 0.0).qv

    ρw = CC.Geometry.WVector.(zeros(FT, face_space))
    ρw0 = 0.0
    term_vel_rai = CC.Geometry.WVector.(zeros(FT, face_space))
    term_vel_sno = CC.Geometry.WVector.(zeros(FT, face_space))

    if moisture isa EquilibriumMoisture
        ts = @. TD.PhaseEquil_ρθq(params, ip.ρ, ip.θ_liq_ice, ip.q_tot)
    elseif moisture isa NonEquilibriumMoisture
        q = @. TD.PhasePartition(ip.q_tot, ip.q_liq, ip.q_ice)
        ts = @. TD.PhaseNonEquil_ρθq(params, ip.ρ, ip.θ_liq_ice, q)
    else
        error(
            "Wrong moisture choise $moisture. The supported options are EquilibriumMoisture and NonEquilibriumMoisture",
        )
    end

    return CC.Fields.FieldVector(;
        constants = (; ρ = ip.ρ, p = ip.p, θ_liq_ice = ip.θ_liq_ice),
        moisture_variables = (;
            T = ip.T,
            θ_dry = ip.θ_dry,
            q_tot = ip.q_tot,
            q_liq = ip.q_liq,
            q_ice = ip.q_ice,
            ts = ts,
        ),
        precip_variables = (; q_rai = ip.q_rai, q_sno = ip.q_sno),
        moisture_sources = (; S_q_liq = ip.S_ql_moisture, S_q_ice = ip.S_qi_moisture),
        precip_sources = (;
            S_q_tot = ip.S_qt_precip,
            S_q_liq = ip.S_ql_precip,
            S_q_ice = ip.S_qi_precip,
            S_q_rai = ip.S_qr_precip,
            S_q_sno = ip.S_qs_precip,
        ),
        precip_velocities = (; term_vel_rai = term_vel_rai, term_vel_sno = term_vel_sno),
        ρw = ρw,
        params = params,
        w_params = w_params,
        q_surf = q_surf,
        ρw0 = ρw0,
        Stats = Stats,
        TS = TS,
    )
end
