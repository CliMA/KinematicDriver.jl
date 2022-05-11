function zero_tendencies!(dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    @. dY.ρq_tot = FT(0)
end

function precompute_aux!(dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    ts = @. TD.PhaseEquil_ρθq(aux.params, aux.ρ, aux.θ_liq_ice, Y.ρq_tot / aux.ρ)

    @. aux.q_liq = TD.liquid_specific_humidity(aux.params, ts)
    @. aux.q_ice = TD.ice_specific_humidity(aux.params, ts)
    @. aux.T     = TD.air_temperature(aux.params, ts)
    @. aux.ρw    = Geometry.WVector.(aux.w_params.w1 * sin(pi * t / aux.w_params.t1))

    aux.ρw0      = aux.w_params.w1 * sin(pi * t / aux.w_params.t1)
end

# Advection Equation: ∂ϕ/dt = -∂(vΦ)
function advection_tendency!(dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    fcc = Operators.FluxCorrectionC2C(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    If = Operators.InterpolateC2F()
    ∂ = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.WVector(aux.ρw0 * aux.q_surf)), # TODO: change to correct ρvw_0 b.c.
        top = Operators.Extrapolate(),
    )
    @. dY.ρq_tot  += -(∂(aux.ρw / If(aux.ρ) * If(Y.ρq_tot))) + fcc(aux.ρw / If(aux.ρ), Y.ρq_tot)

    #Spaces.weighted_dss!(dY.ρq_tot)
    return dY
end

function sources_tendency!(dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    return dY
end

function rhs!(dY, Y, aux, t)

    zero_tendencies!(dY, Y, aux, t)

    precompute_aux!(dY, Y, aux, t)

    advection_tendency!(dY, Y, aux, t)

    sources_tendency!(dY, Y, aux, t)
end
