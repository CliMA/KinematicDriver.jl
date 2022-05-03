function zero_tendencies!(dY, Y, aux, t)
    FT = eltype(Y.q_tot)
    @. dY.q_tot = FT(0)
end

function precompute_aux!(dY, Y, aux, t)
    FT = eltype(Y.q_tot)
    ts = @. TD.PhaseEquil_ρθq(aux.params, aux.ρ, aux.θ_liq_ice, Y.q_tot)

    @. aux.q_liq = TD.liquid_specific_humidity(aux.params, ts)
    @. aux.q_ice = TD.ice_specific_humidity(aux.params, ts)
    @. aux.T     = TD.air_temperature(aux.params, ts)

end

# Advection Equation: ∂ϕ/dt = -∂(vΦ)
function advection_tendency!(dY, Y, aux, t)
    FT = eltype(Y.q_tot)

    # TODO @. w = Y.w * f(t)

    fcc = Operators.FluxCorrectionC2C(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )

    A_qt = Operators.AdvectionC2C(
           bottom = Operators.Extrapolate(),
           top = Operators.Extrapolate(),
    )

<<<<<<< HEAD
    @. dθ = 0.0
    @. dqv = -A_qv(w, qv) + fcc(w, qv)
    return dY
end

#TODO - add microphysics tendency
=======
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
>>>>>>> b5c92ac36415e8fbe849cb749d3e4557ca3bd7d9
