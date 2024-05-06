"""
    Returns the number of new activated aerosol particles and updates aerosol number density
"""
@inline function aerosol_activation_helper(
    kid_params,
    thermo_params,
    air_params,
    activation_params,
    q_tot,
    q_liq,
    N_aer,
    N_aer_0,
    T,
    p,
    ρ,
    ρw,
    dt,
)

    FT = eltype(q_tot)
    S_Nl::FT = FT(0)
    S_Na::FT = FT(0)

    q = TD.PhasePartition(q_tot, q_liq, FT(0))
    S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Liquid())

    if (S < FT(0))
        return (; S_Nl, S_Na)
    end

    (; r_dry, std_dry, κ) = kid_params
    w = ρw / ρ

    aerosol_distribution = CMAM.AerosolDistribution((CMAM.Mode_κ(r_dry, std_dry, N_aer_0, FT(1), FT(1), FT(0), κ, 1),))
    N_act = CMAA.N_activated_per_mode(activation_params, aerosol_distribution, air_params, thermo_params, T, p, w, q)[1]

    if isnan(N_act)
        return (; S_Nl, S_Na)
    end

    S_Nl = max(0, N_act - (N_aer_0 - N_aer)) / dt
    S_Na = -S_Nl

    return (; S_Nl, S_Na)
end

@inline function aerosol_activation_helper(
    kid_params,
    thermo_params,
    air_params,
    activation_params,
    q_tot,
    q_liq,
    N_aer,
    N_aer_0,
    T,
    p,
    ρ,
    ρw,
    dt,
    moments
)
    tmp = aerosol_activation_helper(    
        kid_params,
        thermo_params,
        air_params,
        activation_params,
        q_tot,
        q_liq,
        N_aer,
        N_aer_0,
        T,
        p,
        ρ,
        ρw,
        dt,
    )
    FT = eltype(q_tot)
    S_act = ntuple(length(moments)) do k
        k > 1 ? FT(0) : tmp.S_Nl
    end

    return (; tmp.S_Nl, tmp.S_Na, S_act)
end

"""
Aerosol activation tendencies
"""
@inline function precompute_aux_activation!(sp::CO.AbstractPrecipitationStyle, dY, Y, aux, t)
    error("activation_tendency not implemented for a given $sp")
end
@inline function precompute_aux_activation!(
    ::Union{CO.NoPrecipitation, CO.Precipitation0M, CO.Precipitation1M},
    dY,
    Y,
    aux,
    t,
) end
@inline function precompute_aux_activation!(::CO.Precipitation2M, dY, Y, aux, t)

    aux.aerosol_variables.N_aer = Y.N_aer
    tmp = @. aerosol_activation_helper(
        aux.kid_params,
        aux.thermo_params,
        aux.air_params,
        aux.activation_params,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.aerosol_variables.N_aer,
        aux.aerosol_variables.N_aer_0,
        aux.moisture_variables.T,
        aux.moisture_variables.p,
        aux.moisture_variables.ρ,
        CC.Operators.InterpolateF2C().(aux.prescribed_velocity.ρw.components.data.:1),
        aux.TS.dt,
    )
    aux.activation_sources.S_N_aer = tmp.S_Na
    aux.activation_sources.S_N_liq = tmp.S_Nl
end
@inline function precompute_aux_activation!(::CO.CloudyPrecip, dY, Y, aux, t)

    aux.aerosol_variables.N_aer = Y.N_aer
    tmp = @. aerosol_activation_helper(
        aux.kid_params,
        aux.thermo_params,
        aux.air_params,
        aux.activation_params,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.aerosol_variables.N_aer,
        aux.aerosol_variables.N_aer_0,
        aux.moisture_variables.T,
        aux.moisture_variables.p,
        aux.moisture_variables.ρ,
        CC.Operators.InterpolateF2C().(aux.prescribed_velocity.ρw.components.data.:1),
        aux.TS.dt,
        Y.moments
    )
    aux.activation_sources.S_N_aer = tmp.S_Na
    aux.cloudy_sources.S_activation = tmp.S_act
end

"""
    Prescribed momentum flux as a function of time
"""
@inline function ρw_helper(t, w1, t1)
    return t < t1 ? w1 * sin(pi * t / t1) : 0.0
end

@inline function precompute_aux_prescribed_velocity!(aux, t)

    FT = eltype(aux.moisture_variables.q_tot)
    ρw = FT(ρw_helper(t, aux.kid_params.w1, aux.kid_params.t1))

    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw

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
    @. dY.ρq_tot += -∂(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) * If(Y.ρq_tot))

    if Bool(aux.kid_params.qtot_flux_correction)
        fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
        @. dY.ρq_tot += fcc(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ), Y.ρq_tot)
    end

    return dY
end
@inline function advection_tendency!(::CO.CloudyMoisture, dY, Y, aux, t)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * aux.q_surf)),
        top = CC.Operators.Extrapolate(),
    )
    @. dY.ρq_vap += -∂(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) * If(Y.ρq_vap))

    if Bool(aux.kid_params.qtot_flux_correction)
        fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
        @. dY.ρq_vap += fcc(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ), Y.ρq_vap)
    end

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

    @. dY.ρq_tot += -∂_qt(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) * If(Y.ρq_tot))
    @. dY.ρq_liq += -∂_ql(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) * If(Y.ρq_liq))
    @. dY.ρq_ice += -∂_qi(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) * If(Y.ρq_ice))


    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    if Bool(aux.kid_params.qtot_flux_correction)
        @. dY.ρq_tot += fcc(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ), Y.ρq_tot)
    end
    @. dY.ρq_liq += fcc(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ), Y.ρq_liq)
    @. dY.ρq_ice += fcc(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ), Y.ρq_ice)


    return dY
end
@inline function advection_tendency!(::Union{CO.NoPrecipitation, CO.Precipitation0M}, dY, Y, aux, t) end
@inline function advection_tendency!(::CO.Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(0.0)),
    )

    @. dY.ρq_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
                CC.Geometry.WVector(If(aux.precip_velocities.term_vel_rai) * FT(-1))
            ) * If(Y.ρq_rai),
        )
    @. dY.ρq_sno +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
                CC.Geometry.WVector(If(aux.precip_velocities.term_vel_sno) * FT(-1))
            ) * If(Y.ρq_sno),
        )

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.ρq_rai += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
            CC.Geometry.WVector(If(aux.precip_velocities.term_vel_rai) * FT(-1))
        ),
        Y.ρq_rai,
    )
    @. dY.ρq_sno += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
            CC.Geometry.WVector(If(aux.precip_velocities.term_vel_sno) * FT(-1))
        ),
        Y.ρq_sno,
    )


    return dY
end
@inline function advection_tendency!(::CO.Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(0.0)),
    )

    @. dY.N_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
                CC.Geometry.WVector(If(aux.precip_velocities.term_vel_N_rai) * FT(-1))
            ) * If(Y.N_rai),
        )
    @. dY.ρq_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
                CC.Geometry.WVector(If(aux.precip_velocities.term_vel_rai) * FT(-1))
            ) * If(Y.ρq_rai),
        )

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.N_rai += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
            CC.Geometry.WVector(If(aux.precip_velocities.term_vel_N_rai) * FT(-1))
        ),
        Y.N_rai,
    )
    @. dY.ρq_rai += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
            CC.Geometry.WVector(If(aux.precip_velocities.term_vel_rai) * FT(-1))
        ),
        Y.ρq_rai,
    )

    return dY
end

# TODO: make it work!
@inline function advection_tendency!(::CO.CloudyPrecip, dY, Y, aux, t)
    FT = eltype(Y.ρq_vap)
    Nmom = Int(sum(aux.cloudy_params.NProgMoms))
    
    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(FT(0))),
    )
    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())

    for i in 1:Nmom 
        @. dY.moments.:($$i) = -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
                CC.Geometry.WVector(If(aux.cloudy_velocity.weighted_vt.:($$i)) * FT(-1))
            ) * If(Y.moments.:($$i)),
        )
        @. dY.moments.:($$i) += fcc(
            (
                aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) +
                CC.Geometry.WVector(If(aux.cloudy_velocity.weighted_vt.:($$i)) * FT(-1))
            ),
            Y.moments.:($$i)
        )
    end

    return dY
end
