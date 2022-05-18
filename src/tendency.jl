"""
   Different components of the ODE `rhs!` function,
   depending on the moisture and precipitation types.
"""

"""
    Zero out previous timestep tendencies
"""
function zero_tendencies!(::EquilibriumMoisture, ::NoPrecipitation, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    @. dY.ρq_tot = FT(0)
end

"""
     Precompute the auxiliary values
"""
function precompute_aux!(::EquilibriumMoisture, ::NoPrecipitation, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    ts = @. TD.PhaseEquil_ρθq(aux.params, aux.ρ, aux.θ_liq_ice, Y.ρq_tot / aux.ρ)

    @. aux.q_tot = TD.total_specific_humidity(aux.params, ts)
    @. aux.q_liq = TD.liquid_specific_humidity(aux.params, ts)
    @. aux.q_ice = TD.ice_specific_humidity(aux.params, ts)
    @. aux.T     = TD.air_temperature(aux.params, ts)
    @. aux.ρw    = CC.Geometry.WVector.(aux.w_params.w1 * sin(pi * t / aux.w_params.t1))

    aux.ρw0      = aux.w_params.w1 * sin(pi * t / aux.w_params.t1)
end

"""
   Advection Equation: ∂ϕ/dt = -∂(vΦ)
"""
function advection_tendency!(::EquilibriumMoisture, ::NoPrecipitation, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    fcc = CC.Operators.FluxCorrectionC2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.Extrapolate(),
    )
    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.ρw0 * aux.q_surf)), # TODO: change to correct ρvw_0 b.c.
        top = CC.Operators.Extrapolate(),
    )
    @. dY.ρq_tot  += -(∂(aux.ρw / If(aux.ρ) * If(Y.ρq_tot))) + fcc(aux.ρw / If(aux.ρ), Y.ρq_tot)

    #CC.Spaces.weighted_dss!(dY.ρq_tot)
    return dY
end

"""
   Additional source terms
"""
function sources_tendency!(::EquilibriumMoisture, ::NoPrecipitation, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    @. dY.ρq_tot += FT(0)

    return dY
end
