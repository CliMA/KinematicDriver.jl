function zero_tendencies!(dY, Y, aux, t)
    FT = eltype(Y.q_tot)
    @. dY.q_tot = FT(0)
end

function dρ_dz!(ρ, params_aux_Y, z)

    FT = eltype(ρ)

    params = params_aux_Y.params
    aux = params_aux_Y.aux
    Y = params_aux_Y.Y

    # initial profiles
    θ::FT = aux.θ
    q_tot::FT = Y.q_tot
    q_liq::FT = aux.q_liq
    q_ice::FT = aux.q_ice
    q_v::FT = q_tot - q_liq - q_ice

    q = TD.PhasePartition(q_v, q_liq, q_ice)

    # constants
    g::FT = CP.Planet.grav(params)
    cp_m::FT = TD.cp_m(params, q)
    R_m::FT = TD.gas_constant_air(params, q)

    T::FT = θ * (ρ * θ / CP.Planet.MSLP(params) * R_m)^((R_m / cp_m) / (1 - R_m / cp_m))

    return g / T * ρ * (R_m / cp_m - 1) / R_m
end

function ρ_ivp(::Type{FT}, params, aux, Y) where {FT}

    init_surface = init_condition(FT, params, 0.0)

    ρ_0::FT = init_surface.ρ_0
    z_0::FT = init_surface.z_0
    z_max::FT = init_surface.z_2

    z_span = (z_0, z_max)
    params_aux_Y = (params=params, aux=aux, Y=Y)
    prob = ODEProblem(dρ_dz!, ρ_0, z_span, params_aux_Y)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    return sol
end

function precompute_aux!(dY, Y, aux, t)
    FT = eltype(Y.q_tot)
    # TODO: update density
    ρ_profile = ρ_ivp(FT, params, aux, Y)
    @. aux.ρ = ρ_profile(aux.coord.z)

    ts = @. TD.PhaseEquil_ρθq(aux.params, aux.ρ, aux.θ_liq_ice, Y.q_tot)

    @. aux.q_liq = TD.liquid_specific_humidity(aux.params, ts)
    @. aux.q_ice = TD.ice_specific_humidity(aux.params, ts)
    @. aux.T     = TD.air_temperature(aux.params, ts)
    @. aux.ρw     = Geometry.WVector.(aux.w_params.w1 * sin(pi * t / aux.w_params.t1))
end

# Advection Equation: ∂ϕ/dt = -∂(vΦ)
function advection_tendency!(dY, Y, aux, t)
    FT = eltype(Y.q_tot)

    fcc = Operators.FluxCorrectionC2C(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )

    A_qt = Operators.AdvectionC2C(
           bottom = Operators.Extrapolate(),
           top = Operators.Extrapolate(),
    )
    #aux.w = aux.ρw ./ aux.ρ
    @. dY.q_tot = -A_qt(aux.w, Y.q_tot) + fcc(aux.w, Y.q_tot)
    return dY
end

function sources_tendency!(dY, Y, aux, t)
    FT = eltype(Y.q_tot)

    return dY
end

function rhs!(dY, Y, aux, t)

    zero_tendencies!(dY, Y, aux, t)

    precompute_aux!(dY, Y, aux, t)

    advection_tendency!(dY, Y, aux, t)

    sources_tendency!(dY, Y, aux, t)
end
