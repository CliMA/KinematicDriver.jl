"""
  Main building blocks of the KiD model
"""

"""
   Interface to ClimaCore.jl Returns an instance of space and face space
   (which are a discretized function space over our computational domain).
"""
function make_function_space(z_min, z_max, n_elem)
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
function make_rhs_function(
    moisture::AbstractMoistureStyle,
    precipitation::AbstractPrecipitationStyle,
)
    function rhs!(dY, Y, aux, t)
        zero_tendencies!(moisture, precipitation, dY, Y, aux, t)
        precompute_aux!(moisture, precipitation, dY, Y, aux, t)
        advection_tendency!(moisture, precipitation, dY, Y, aux, t)
        sources_tendency!(moisture, precipitation, dY, Y, aux, t)
    end

    return rhs!
end

"""
    Interface to ODE solver. Returns the ClimaCore.jl FieldVector with
    ODE solver state variables. The state is created via dispatching
    on different moisture and precipitation types
"""
function initialise_state(
    ::EquilibriumMoisture,
    ::Union{NoPrecipitation, Precipitation0M},
    initial_profiles,
)
    return CC.Fields.FieldVector(; ρq_tot = initial_profiles.ρq_tot)
end
function initialise_state(
    ::NonEquilibriumMoisture,
    ::Union{NoPrecipitation, Precipitation0M},
    initial_profiles,
)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_liq = initial_profiles.ρq_liq,
        ρq_ice = initial_profiles.ρq_ice,
    )
end
function initialise_state(
    ::EquilibriumMoisture,
    ::Precipitation1M,
    initial_profiles,
)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_rai = initial_profiles.ρq_rai,
        ρq_sno = initial_profiles.ρq_sno,
    )
end
function initialise_state(
    ::NonEquilibriumMoisture,
    ::Precipitation1M,
    initial_profiles,
)
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
function initialise_aux(
    initial_profiles,
    params,
    w_params,
    q_surf,
    ρw0,
    TS,
    Stats,
    face_space,
)
    ρw = CC.Geometry.WVector.(zeros(FT, face_space))

    return CC.Fields.FieldVector(;
        ρ = initial_profiles.ρ,
        T = initial_profiles.T,
        p = initial_profiles.p,
        θ_liq_ice = initial_profiles.θ_liq_ice,
        θ_dry = initial_profiles.θ_dry,
        q_tot = initial_profiles.q_tot,
        q_liq = initial_profiles.q_liq,
        q_ice = initial_profiles.q_ice,
        q_rai = initial_profiles.q_rai,
        q_sno = initial_profiles.q_sno,
        ρw = ρw,
        params = params,
        w_params = w_params,
        q_surf = q_surf,
        ρw0 = ρw0,
        Stats = Stats,
        TS = TS,
    )
end
