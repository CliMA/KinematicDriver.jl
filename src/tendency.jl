# Advection Equation: ∂ϕ/dt = -∂(vΦ)
function advection_tendency!(dY, Y, _, t)

    Yc = Y.Yc
    w = Y.w
    # TODO @. w = Y.w * sin(t) ?

    dYc = dY.Yc

    θ = Yc.θ
    qv = Yc.qv

    dθ = dYc.θ
    dqv = dYc.qv

    fcc = Operators.FluxCorrectionC2C(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    fcf = Operators.FluxCorrectionF2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    A_θ = Operators.AdvectionC2C(
          bottom = Operators.SetValue(279.9),   #TODO - boundary conditions!
          top = Operators.Extrapolate(),
    )

    A_qv = Operators.AdvectionC2C(
           bottom = Operators.SetValue(0.016),   #TODO - boundary conditions!
           top = Operators.Extrapolate(),
    )

    @. dθ = -A_θ(w, θ) + fcc(w, θ)
    @. dqv = -A_qv(w, qv) + fcc(w, qv)
    return dY
end

#TODO - add microphysics tendency
