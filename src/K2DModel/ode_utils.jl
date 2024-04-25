"""
  Main buliding blocks for the ODE solver
"""

"""
   Interface to ClimaCore.jl Returns an instance of space and face space
   (which are a discretized function space over our computational domain).
"""
function make_function_space(FT; xlim = (0.0, 3000.0), zlim = (0.0, 3000.0), helem = 32, velem = 32, npoly = 4)
    vertdomain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(zlim[1]),
        CC.Geometry.ZPoint{FT}(zlim[2]);
        boundary_names = (:bottom, :top),
    )
    vertmesh = CC.Meshes.IntervalMesh(vertdomain, nelems = velem)
    vert_center_space = CC.Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain =
        CC.Domains.IntervalDomain(CC.Geometry.XPoint{FT}(xlim[1]), CC.Geometry.XPoint{FT}(xlim[2]), periodic = true)
    horzmesh = CC.Meshes.IntervalMesh(horzdomain; nelems = helem)
    horztopology = CC.Topologies.IntervalTopology(horzmesh)

    quad = CC.Spaces.Quadratures.GLL{npoly + 1}()
    horzspace = CC.Spaces.SpectralElementSpace1D(horztopology, quad)

    hv_center_space = CC.Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
    hv_face_space = CC.Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)
    return (hv_center_space, hv_face_space)
end

"""
   Interface to ODE solver. Returns the function needed to compute the
   right hand side of the solved ODE. The rhs is assembled via dispatch
   based on the moisture and precipitation types.
"""
function make_rhs_function(ms::CO.AbstractMoistureStyle, ps::CO.AbstractPrecipitationStyle)
    function rhs!(dY, Y, aux, t)

        CO.zero_tendencies!(dY)

        precompute_aux_prescribed_velocity!(aux, t)
        CO.precompute_aux_thermo!(ms, Y, aux)
        K1D.precompute_aux_activation!(ps, dY, Y, aux, t)
        CO.precompute_aux_precip!(ps, Y, aux)

        CO.cloud_sources_tendency!(ms, dY, Y, aux, t)
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
    domain_width,
    domain_height,
    TS,
    Stats,
    space,
    face_space,
    moisture,
    precip,
)

    q_surf = CO.init_profile(FT, kid_params, thermo_params, 0.0).qv

    ρu = CC.Geometry.UVector.(zeros(FT, space))
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
        ),
        (;
            prescribed_velocity = CC.Fields.FieldVector(; ρu = ρu, ρw = ρw, ρw0 = ρw0),
            kid_params = kid_params,
            domain_width = domain_width,
            domain_height = domain_height,
            q_surf = q_surf,
        ),
    )
end
