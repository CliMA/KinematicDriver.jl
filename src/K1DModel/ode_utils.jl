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
        boundary_tags = (:bottom, :top),
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

        for eq_style in [ms, ps]
            CO.zero_tendencies!(eq_style, dY, Y, aux, t)
        end

        precompute_aux_prescribed_velocity!(aux, t)
        CO.precompute_aux_thermo!(ms, dY, Y, aux, t)
        CO.precompute_aux_moisture_sources!(ms, dY, Y, aux, t)
        precompute_aux_activation!(ps, dY, Y, aux, t)
        CO.precompute_aux_precip!(ps, dY, Y, aux, t)

        for eq_style in [ms, ps]
            advection_tendency!(eq_style, dY, Y, aux, t)
            CO.sources_tendency!(eq_style, dY, Y, aux, t)
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

        for eq_style in [ms, ps]
            CO.zero_tendencies!(eq_style, dY, Y, aux, t)
        end

        CO.precompute_aux_thermo!(ms, dY, Y, aux, t)
        CO.precompute_aux_precip!(ps, dY, Y, aux, t)

        for eq_style in [ms, ps]
            advection_tendency!(eq_style, dY, Y, aux, t)
            CO.sources_tendency!(eq_style, dY, Y, aux, t)
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
)

    q_surf = CO.init_profile(FT, kid_params, thermo_params, 0.0).qv

    ρw = CC.Geometry.WVector.(zeros(FT, face_space))
    ρw0 = 0.0

    return merge(
        CO.initialise_aux(FT, ip, common_params, thermo_params, air_params, activation_params, TS, Stats, moisture), 
        (;
            prescribed_velocity = CC.Fields.FieldVector(; ρw = ρw, ρw0 = ρw0),
            kid_params = kid_params,
            q_surf = q_surf,
        )
    )
end
