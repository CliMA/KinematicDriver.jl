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

        for eq_style in [ms, ps]
            K1D.zero_tendencies!(eq_style, dY, Y, aux, t)
        end

        precompute_aux_prescribed_velocity!(aux, t)
        K1D.precompute_aux_thermo!(ms, dY, Y, aux, t)
        K1D.precompute_aux_precip!(ps, dY, Y, aux, t)

        for eq_style in [ms, ps]
            advection_tendency!(eq_style, dY, Y, aux, t)
            K1D.sources_tendency!(eq_style, dY, Y, aux, t)
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
    ic,
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
)

    q_surf = K1D.init_condition(FT, kid_params, thermo_params, 0.0).qv

    ρu = CC.Geometry.UVector.(zeros(FT, space))
    ρw = CC.Geometry.WVector.(zeros(FT, face_space))
    ρw0 = 0.0
    term_vel_rai = CC.Geometry.WVector.(zeros(FT, face_space))
    term_vel_sno = CC.Geometry.WVector.(zeros(FT, face_space))

    if moisture isa CO.EquilibriumMoisture
        ts = @. TD.PhaseEquil_ρθq(thermo_params, ic.ρ, ic.θ_liq_ice, ic.q_tot)
    elseif moisture isa CO.NonEquilibriumMoisture
        q = @. TD.PhasePartition(ic.q_tot, ic.q_liq, ic.q_ice)
        ts = @. TD.PhaseNonEquil_ρθq(thermo_params, ic.ρ, ic.θ_liq_ice, q)
    else
        error(
            "Wrong moisture choise $moisture. The supported options are EquilibriumMoisture and NonEquilibriumMoisture",
        )
    end

    return (;
        moisture_variables = CC.Fields.FieldVector(;
            ρ = ic.ρ,
            ρ_dry = ic.ρ_dry,
            p = ic.p,
            T = ic.T,
            θ_liq_ice = ic.θ_liq_ice,
            θ_dry = ic.θ_dry,
            q_tot = ic.q_tot,
            q_liq = ic.q_liq,
            q_ice = ic.q_ice,
            ts = ts,
        ),
        precip_variables = CC.Fields.FieldVector(;
            q_rai = ic.q_rai,
            q_sno = ic.q_sno,
            N_liq = ic.N_liq,
            N_rai = ic.N_rai,
        ),
        aerosol_variables = CC.Fields.FieldVector(; N_aer = ic.N_aer, N_aer_0 = ic.N_aer_0, S_N_aer = ic.S_Na),
        moisture_sources = CC.Fields.FieldVector(; S_q_liq = ic.S_ql_moisture, S_q_ice = ic.S_qi_moisture),
        precip_sources = CC.Fields.FieldVector(;
            S_q_tot = ic.S_qt_precip,
            S_q_liq = ic.S_ql_precip,
            S_q_ice = ic.S_qi_precip,
            S_q_rai = ic.S_qr_precip,
            S_q_sno = ic.S_qs_precip,
            S_N_liq = ic.S_Nl_precip,
            S_N_rai = ic.S_Nr_precip,
        ),
        precip_velocities = CC.Fields.FieldVector(;
            term_vel_rai = term_vel_rai,
            term_vel_sno = term_vel_sno,
            term_vel_N_rai = term_vel_rai,
        ),
        prescribed_velocity = CC.Fields.FieldVector(; ρu = ρu, ρw = ρw, ρw0 = ρw0),
        thermo_params = thermo_params,
        kid_params = kid_params,
        air_params = air_params,
        activation_params = activation_params,
        domain_width = domain_width,
        domain_height = domain_height,
        q_surf = q_surf,
        Stats = Stats,
        TS = TS,
    )
end
