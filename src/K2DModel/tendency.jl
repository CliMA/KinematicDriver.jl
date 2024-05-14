"""
    Prescribed momentum flux as a function of time
"""
@inline function prescribed_momentum_helper(x, z, t, w1, t1, width, height)
    time_scale = t < t1 ? sin(pi * t / t1) : 0.0
    ρu = 0.5 * w1 * width / height * cos(pi / height * z) * cos(2 * pi / width * x) * time_scale
    ρw = w1 * sin(pi / height * z) * sin(2 * pi / width * x) * time_scale
    return (; ρu, ρw)
end
@inline function precompute_aux_prescribed_velocity!(aux, t)

    FT = eltype(aux.thermo_params)
    coords = CC.Fields.coordinate_field(axes(aux.prescribed_velocity.ρu))
    face_coords = CC.Fields.coordinate_field(axes(aux.prescribed_velocity.ρw))
    aux.prescribed_velocity.ρu = map(
        coord -> CC.Geometry.UVector(
            prescribed_momentum_helper(
                coord.x,
                coord.z,
                t,
                aux.kid_params.w1,
                aux.kid_params.t1,
                aux.domain_width,
                aux.domain_height,
            ).ρu,
        ),
        coords,
    )
    aux.prescribed_velocity.ρw = map(
        coord -> CC.Geometry.WVector(
            prescribed_momentum_helper(
                coord.x,
                coord.z,
                t,
                aux.kid_params.w1,
                aux.kid_params.t1,
                aux.domain_width,
                aux.domain_height,
            ).ρw,
        ),
        face_coords,
    )

    aux.prescribed_velocity.ρw0 = FT(0)

end

"""
   Advection Equation: ∂ϕ/dt = -∂(vΦ)
"""
@inline function advection_tendency!(sm::CO.AbstractMoistureStyle, dY, Y, aux, t)
    error("advection_tendency not implemented for a given $sm")
end
@inline function advection_tendency!(sp::CO.AbstractPrecipitationStyle, dY, Y, aux, t)
    error("advection_tendency not implemented for a given $sp")
end
@inline function advection_tendency!(::CO.EquilibriumMoisture, dY, Y, aux, t)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * aux.q_surf)),
        top = CC.Operators.Extrapolate(),
    )
    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    hdiv = CC.Operators.WeakDivergence()

    @. dY.ρq_tot += -∂(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_tot))
    @. dY.ρq_tot += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_tot)
    @. dY.ρq_tot += -hdiv(aux.prescribed_velocity.ρu / aux.thermo_variables.ρ * Y.ρq_tot)
    CC.Spaces.weighted_dss!(dY.ρq_tot)

    return dY
end
@inline function advection_tendency!(::CO.NonEquilibriumMoisture, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()

    ∂_qt = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * aux.q_surf)),
        top = CC.Operators.Extrapolate(),
    )
    ∂_ql = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * FT(0.0))),
        top = CC.Operators.Extrapolate(),
    )
    ∂_qi = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * FT(0.0))),
        top = CC.Operators.Extrapolate(),
    )

    @. dY.ρq_tot += -∂_qt(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_tot))
    @. dY.ρq_liq += -∂_ql(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_liq))
    @. dY.ρq_ice += -∂_qi(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_ice))


    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.ρq_tot += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_tot)
    @. dY.ρq_liq += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_liq)
    @. dY.ρq_ice += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_ice)

    hdiv = CC.Operators.WeakDivergence()
    @. dY.ρq_tot += -hdiv(aux.prescribed_velocity.ρu / aux.thermo_variables.ρ * Y.ρq_tot)
    @. dY.ρq_liq += -hdiv(aux.prescribed_velocity.ρu / aux.thermo_variables.ρ * Y.ρq_liq)
    @. dY.ρq_ice += -hdiv(aux.prescribed_velocity.ρu / aux.thermo_variables.ρ * Y.ρq_ice)
    CC.Spaces.weighted_dss!(dY.ρq_tot)
    CC.Spaces.weighted_dss!(dY.ρq_liq)
    CC.Spaces.weighted_dss!(dY.ρq_ice)


    return dY
end
@inline function advection_tendency!(::Union{CO.NoPrecipitation, CO.Precipitation0M}, dY, Y, aux, t) end
@inline function advection_tendency!(::CO.Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())

    @. dY.ρq_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
            ) * If(Y.ρq_rai),
        )
    @. dY.ρq_sno +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_sno) * FT(-1))
            ) * If(Y.ρq_sno),
        )

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.ρq_rai += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
        ),
        Y.ρq_rai,
    )
    @. dY.ρq_sno += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_sno) * FT(-1))
        ),
        Y.ρq_sno,
    )

    hdiv = CC.Operators.WeakDivergence()
    @. dY.ρq_rai += -hdiv(aux.prescribed_velocity.ρu / aux.thermo_variables.ρ * Y.ρq_rai)
    @. dY.ρq_sno += -hdiv(aux.prescribed_velocity.ρu / aux.thermo_variables.ρ * Y.ρq_sno)
    CC.Spaces.weighted_dss!(dY.ρq_rai)
    CC.Spaces.weighted_dss!(dY.ρq_sno)

    return dY
end
@inline function advection_tendency!(::CO.Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())

    @. dY.N_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_N_rai) * FT(-1))
            ) * If(Y.N_rai),
        )
    @. dY.ρq_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
            ) * If(Y.ρq_rai),
        )

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.N_rai += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_N_rai) * FT(-1))
        ),
        Y.N_rai,
    )
    @. dY.ρq_rai += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
        ),
        Y.ρq_rai,
    )


    hdiv = CC.Operators.WeakDivergence()
    @. dY.N_rai += -hdiv(aux.prescribed_velocity.ρu / aux.thermo_variables.ρ * Y.N_rai)
    @. dY.ρq_rai += -hdiv(aux.prescribed_velocity.ρu / aux.thermo_variables.ρ * Y.ρq_rai)
    CC.Spaces.weighted_dss!(dY.ρq_rai)
    CC.Spaces.weighted_dss!(dY.N_rai)

    return dY
end
