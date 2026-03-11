"""
  Main buliding blocks for the ODE solver
"""

"""
   Interface to ClimaCore.jl Returns an instance of space and face space
   (which are a discretized function space over our computational domain).
"""
function make_function_space(FT, z_min, z_max, n_elem)
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(z_min),
        CC.Geometry.ZPoint{FT}(z_max),
        boundary_names = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems = n_elem)

    space = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    face_space = CC.Spaces.FaceFiniteDifferenceSpace(space)

    return space, face_space
end

"""
   Interface to ODE solver. Returns the function needed to compute the
   right hand side of the solved ODE for advection, condensation,
   collision, sedimentation and evaporation processes. The rhs is
   assembled via dispatch based on the moisture and precipitation types.
"""
function make_rhs_function(ms::CO.AbstractMoistureStyle, ps::CO.AbstractPrecipitationStyle)
    function rhs!(dY, Y, aux, t)
        CO.zero_tendencies!(dY)

        # defined in `K1DModel/tendency.jl`
        # calc: aux.prescribed_velocity.ρw
        precompute_aux_prescribed_velocity!(aux, t)
        # defined in `Common/tendency.jl`
        # For NonEquilibriumMoisture:
        # calc: aux.thermo_variables.{ts, ρ, p, θ_dry, θ_liq_ice}
        # calc: aux.microph_variables.{q_tot, q_liq, q_ice}
        CO.precompute_aux_thermo!(ms, ps, Y, aux)
        # defined in `Common/tendency.jl`
        # For Precipitation2M:
        # calc: aux.microph_variables.{q_rai, q_sno, N_rai, N_liq}  <-- NOTE: You can assume `q_sno = ρq_sno = 0`
        # calc: aux.velocities.{term_vel_rai, term_vel_N_rai}
        # For IcePrecipitationP3:
        # calc: aux.microph_variables.{ρq_ice, ρq_rim, N_ice, B_rim} <-- maybe not needed
        # calc: aux.velocities.{term_vel_ice, term_vel_N_ice}
        CO.precompute_aux_precip!(ps, Y, aux)  # <-- calls fn for Precipitation2M

        # defined in `K1DModel/tendency.jl`
        # calc: aux.activation_sources.{N_aer, N_liq}
        precompute_aux_activation!(ps, dY, Y, aux, t)  # <-- calls fn for Precipitation2M

        # defined in `Common/tendency.jl`
        # calls: `precompute_aux_moisture_sources` (defined in `Common/tendency.jl`)
        #   calc: aux.cloud_sources.{q_liq, q_ice}
        # increments: dY.{ρq_liq, ρq_ice} += ρ * aux.cloud_sources.{q_...}
        CO.cloud_sources_tendency!(ms, ps, dY, Y, aux, t)
        # defined in `Common/tendency.jl`
        # For Precipitation2M:
        # calls: `precompute_aux_precip_sources` (defined in `Common/tendency.jl`)
        #   calc: aux.precip_sources.{q_tot, q_liq, q_rai, N_liq, N_rai, N_aer}
        # increments: dY.{ρq_tot, ρq_rai, ρq_liq, N_liq, N_rai, N_aer}
        #               += aux.precip_sources.{...} 
        #        and/or += aux.activation_sources.{N_liq, N_aer}
        # For Precipitation2M_P3:
        # calls: `precompute_aux_precip_sources` (defined in `Common/tendency.jl`)
        # TODO: write collisions
        CO.precip_sources_tendency!(ms, ps, dY, Y, aux, t)  # <-- calls fn for Precipitation2M + write code

        for eq_style in [ms, ps]
            # defined in  `K1DModel/tendency.jl`
            # For Precipitation2M:
            # increments: dY.{ρq_rai, N_liq, N_rai, N_aer} += ... * aux.velocities.{term_vel_N_rai, term_vel_rai}
            # For NonEquilibriumMoisture:
            # increments: dY.{ρq_tot, ρq_liq, ρq_ice} += {aux.prescribed_velocity} * {Y.ρq_...}
            # For Precipitation2M_P3
            # - delegates to `Precipitation2M` (see above) and `IcePrecipitationP3`
            # For IcePrecipitationP3:
            # increments: dY.{N_ice, B_rim, ρq_rim} += ... * aux.velocities.{term_vel_N_ice, term_vel_q_ice}
            advection_tendency!(eq_style, dY, Y, aux, t)
        end

    end
    return rhs!
end

"""
   Interface to ODE solver. Returns the function needed to compute the
   right hand side of the solved ODE for collision and sedimentation
   processes. The rhs is assembled via dispatch based on precipitation type.
"""
function make_rhs_function_col_sed(ms::CO.AbstractMoistureStyle, ps::CO.AbstractPrecipitationStyle)

    function rhs!(dY, Y, aux, t)

        CO.zero_tendencies!(dY)

        CO.precompute_aux_thermo!(ms, ps, Y, aux)
        CO.precompute_aux_precip!(ps, Y, aux)

        CO.precip_sources_tendency!(ms, ps, dY, Y, aux, t)

        for eq_style in [ms, ps]
            advection_tendency!(eq_style, dY, Y, aux, t)
        end

    end
    return rhs!
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
    kid_params,
    thermo_params,
    air_params,
    activation_params,
    TS,
    Stats,
    face_space,
    moisture,
    precip,
    cloudy_params = nothing,
)
    q_surf = CO.init_profile(FT, kid_params, thermo_params, 0.0).qv

    ρw = CC.Geometry.WVector.(zeros(FT, face_space))
    ρw0 = 0.0

    return merge(
        CO.initialise_aux(
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
            cloudy_params,
        ),
        (; prescribed_velocity = CC.Fields.FieldVector(; ρw = ρw, ρw0 = ρw0), kid_params = kid_params, q_surf = q_surf),
    )
end
