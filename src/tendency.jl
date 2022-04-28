# Advection Equation: ∂ϕ/dt = -∂(vΦ)
function advection_tendency!(dY, Y, _, t, w_params)

    Yc = Y.Yc
    w = Y.w
    w1 = w_params.w1
    t1 = w_params.t1
    w = t < t1 ? w1 * sin.(π * t / t1) : 0
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

    @. dθ = 0.0
    @. dqv = -A_qv(w, qv) + fcc(w, qv)
    return dY
end

#TODO - add microphysics tendency