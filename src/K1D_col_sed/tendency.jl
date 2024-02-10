"""
   Different components of the ODE `rhs!` function,
   depending on the moisture and precipitation types.
"""

"""
    Zero out previous timestep tendencies and aux.S source terms
"""
@inline function zero_tendencies!(sp::K1D.AbstractPrecipitationStyle, dY, Y, aux, t)
    error("zero_tendecies not implemented for a given $sp")
end
@inline function zero_tendencies!(::Union{K1D.NoPrecipitation, K1D.Precipitation0M}, dY, Y, aux, t) end
@inline function zero_tendencies!(::K1D.Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    @. dY.ρq_tot = FT(0)
    @. dY.ρq_liq = FT(0)
    @. dY.ρq_rai = FT(0)
end
@inline function zero_tendencies!(::K1D.Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    @. dY.ρq_tot = FT(0)
    @. dY.ρq_liq = FT(0)
    @. dY.ρq_rai = FT(0)
    @. dY.N_liq = FT(0)
    @. dY.N_rai = FT(0)
end

"""
    helper functions
"""
@inline function update_density!(Y, aux)
    # here we assume vapor is already considered in rhod
    aux.moisture_variables.ρ = aux.moisture_variables.ρ_dry .+ Y.ρq_tot
end

@inline function precip_helper_sources!(ps::K1D.Precipitation0M, q_liq, dt)

    q = TD.PhasePartition(q_liq, q_liq, 0.0)
    S_qt = -min(max(0, (q_liq) / dt), -CM0.remove_precipitation(ps.params, q))
    S_q_rai = -S_qt
    S_q_liq = S_qt
    S_q_tot = S_qt

    return (; S_q_tot, S_q_liq, S_q_rai)
end
@inline function precip_helper_sources!(ps::K1D.Precipitation1M, kid_params, q_liq, q_rai, ρ, dt)

    FT = eltype(q_liq)
    S_q_tot::FT = FT(0)
    S_q_liq::FT = FT(0)
    S_q_rai::FT = FT(0)

    if Bool(kid_params.precip_sources)
        tmp::FT = FT(0)
        rf = ps.rain_formation
        if rf isa CMP.Acnv1M{FT}
            tmp = CM1.conv_q_liq_to_q_rai(rf, q_liq, smooth_transition = true)
        elseif typeof(rf) in [CMP.LD2004{FT}, CMP.VarTimescaleAcnv{FT}]
            tmp = CM2.conv_q_liq_to_q_rai(rf, q_liq, ρ; N_d = kid_params.prescribed_Nd)
        elseif typeof(rf.acnv) in [CMP.AcnvKK2000{FT}, CMP.AcnvB1994{FT}, CMP.AcnvTC1980{FT}]
            tmp = CM2.conv_q_liq_to_q_rai(rf, q_liq, ρ; N_d = kid_params.prescribed_Nd)
        else
            error("Unrecognized rain formation scheme")
        end
        S_qt_rain = -min(max(FT(0), q_liq / dt), tmp)
        S_q_tot += S_qt_rain
        S_q_liq += S_qt_rain
        S_q_rai -= S_qt_rain

        # accretion cloud water + rain
        if typeof(rf) in [CMP.Acnv1M{FT}, CMP.LD2004{FT}, CMP.VarTimescaleAcnv{FT}]
            tmp = CM1.accretion(ps.liquid, ps.rain, ps.sedimentation.rain, ps.ce, q_liq, q_rai, ρ)
        elseif typeof(rf.accr) in [CMP.AccrKK2000{FT}, CMP.AccrB1994{FT}]
            tmp = CM2.accretion(rf, q_liq, q_rai, ρ)
        elseif rf.accr isa CMP.AccrTC1980{FT}
            tmp = CM2.accretion(rf, q_liq, q_rai)
        else
            error("Unrecognized rain formation scheme")
        end
        S_qr = min(max(FT(0), q_liq / dt), tmp)
        S_q_tot -= S_qr
        S_q_liq -= S_qr
        S_q_rai += S_qr
    end

    term_vel_rai = CM1.terminal_velocity(ps.rain, ps.sedimentation.rain, ρ, q_rai)

    return (; S_q_tot, S_q_liq, S_q_rai, term_vel_rai)
end
@inline function precip_helper_sources!(ps::K1D.Precipitation2M, kid_params, q_liq, q_rai, N_liq, N_rai, ρ, dt)

    sb2006 = ps.rain_formation
    vel_scheme = ps.sedimentation

    FT = eltype(q_liq)

    S_q_tot::FT = FT(0)
    S_q_liq::FT = FT(0)
    S_q_rai::FT = FT(0)
    S_N_liq::FT = FT(0)
    S_N_rai::FT = FT(0)

    if Bool(kid_params.precip_sources)

        # autoconversion liquid to rain
        tmp = CM2.autoconversion(sb2006.acnv, q_liq, q_rai, ρ, N_liq)
        S_qr = min(max(FT(0), q_liq / dt), tmp.dq_rai_dt)
        S_q_tot -= S_qr
        S_q_liq -= S_qr
        S_q_rai += S_qr
        S_Nr = min(max(FT(0), N_liq / 2 / dt), tmp.dN_rai_dt)
        S_N_liq -= 2 * S_Nr
        S_N_rai += S_Nr

        # self_collection
        tmp_l = CM2.liquid_self_collection(sb2006.acnv, q_liq, ρ, -2 * S_Nr)
        tmp_r = CM2.rain_self_collection(sb2006.pdf, sb2006.self, q_rai, ρ, N_rai)
        S_N_liq += -min(max(FT(0), N_liq / dt), -tmp_l)
        S_N_rai += -min(max(FT(0), N_rai / dt), -tmp_r)

        # accretion cloud water + rain
        tmp = CM2.accretion(sb2006, q_liq, q_rai, ρ, N_liq)
        S_qr = min(max(FT(0), q_liq / dt), tmp.dq_rai_dt)
        S_q_tot -= S_qr
        S_q_liq -= S_qr
        S_q_rai += S_qr
        S_N_liq += -min(max(FT(0), N_liq / dt), -tmp.dN_liq_dt)
        S_N_rai += FT(0)
    end

    if Bool(kid_params.precip_sinks)

        # rain breakup
        tmp = CM2.rain_breakup(sb2006.pdf, sb2006.brek, q_rai, ρ, N_rai, tmp_r)
        S_N_rai += min(max(FT(0), N_rai / dt), tmp_r)

    end

    vt = CM2.rain_terminal_velocity(sb2006, vel_scheme, q_rai, ρ, N_rai)
    term_vel_N_rai = vt[1]
    term_vel_rai = vt[2]

    return (; S_q_tot, S_q_liq, S_q_rai, S_N_liq, S_N_rai, term_vel_rai, term_vel_N_rai)
end

"""
     Precompute the auxiliary values
"""
@inline function precompute_aux_precip!(sp::K1D.AbstractPrecipitationStyle, dY, Y, aux, t)
    error("precompute_aux not implemented for a given $sp")
end
@inline function precompute_aux_precip!(::K1D.NoPrecipitation, dY, Y, aux, t)
    FT = eltype(Y.ρq_liq)

    aux.precip_variables.q_rai = FT(0)

    aux.precip_sources.S_q_tot = FT(0)
    aux.precip_sources.S_q_liq = FT(0)
    aux.precip_sources.S_q_rai = FT(0)
end
@inline function precompute_aux_precip!(ps::K1D.Precipitation0M, dY, Y, aux, t)
    FT = eltype(Y.ρq_liq)

    aux.precip_variables.q_rai = FT(0)

    tmp = @. precip_helper_sources!(ps, aux.moisture_variables.q_liq, aux.TS.dt)
    aux.precip_sources.S_q_tot = tmp.S_q_tot
    aux.precip_sources.S_q_liq = tmp.S_q_liq
    aux.precip_sources.S_q_rai = tmp.S_q_rai
end
@inline function precompute_aux_precip!(ps::K1D.Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_liq)

    aux.moisture_variables.q_liq = Y.ρq_liq ./ aux.moisture_variables.ρ
    aux.precip_variables.q_rai = Y.ρq_rai ./ aux.moisture_variables.ρ

    tmp = @. precip_helper_sources!(
        ps,
        aux.kid_params,
        aux.moisture_variables.q_liq,
        aux.precip_variables.q_rai,
        aux.moisture_variables.ρ,
        aux.TS.dt,
    )
    aux.precip_sources.S_q_tot = tmp.S_q_tot
    aux.precip_sources.S_q_liq = tmp.S_q_liq
    aux.precip_sources.S_q_rai = tmp.S_q_rai

    If = CC.Operators.InterpolateC2F(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. aux.precip_velocities.term_vel_rai = CC.Geometry.WVector(If(tmp.term_vel_rai) * FT(-1))
end
@inline function precompute_aux_precip!(ps::K1D.Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_liq)

    aux.moisture_variables.q_liq = Y.ρq_liq ./ aux.moisture_variables.ρ
    aux.precip_variables.q_rai = Y.ρq_rai ./ aux.moisture_variables.ρ
    aux.precip_variables.N_liq = Y.N_liq
    aux.precip_variables.N_rai = Y.N_rai
    tmp = @. precip_helper_sources!(
        ps,
        aux.kid_params,
        aux.moisture_variables.q_liq,
        aux.precip_variables.q_rai,
        aux.precip_variables.N_liq,
        aux.precip_variables.N_rai,
        aux.moisture_variables.ρ,
        aux.TS.dt,
    )
    aux.precip_sources.S_q_tot = tmp.S_q_tot
    aux.precip_sources.S_q_liq = tmp.S_q_liq
    aux.precip_sources.S_q_rai = tmp.S_q_rai
    aux.precip_sources.S_N_liq = tmp.S_N_liq
    aux.precip_sources.S_N_rai = tmp.S_N_rai

    If = CC.Operators.InterpolateC2F(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. aux.precip_velocities.term_vel_N_rai = CC.Geometry.WVector(If(tmp.term_vel_N_rai) * FT(-1))
    @. aux.precip_velocities.term_vel_rai = CC.Geometry.WVector(If(tmp.term_vel_rai) * FT(-1))
end

"""
   Advection Equation: ∂ϕ/dt = -∂(vΦ)
"""
@inline function advection_tendency!(sp::K1D.AbstractPrecipitationStyle, dY, Y, aux, t)
    error("advection_tendency not implemented for a given $sp")
end
@inline function advection_tendency!(::Union{K1D.NoPrecipitation, K1D.Precipitation0M}, dY, Y, aux, t) end
@inline function advection_tendency!(::K1D.Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_liq)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(0.0)),
    )
    @. dY.ρq_rai += -∂(aux.precip_velocities.term_vel_rai * If(Y.ρq_rai))

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.ρq_rai += fcc(aux.precip_velocities.term_vel_rai, Y.ρq_rai)

    return dY
end
@inline function advection_tendency!(::K1D.Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_liq)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(0.0)),
    )
    @. dY.N_rai += -∂(aux.precip_velocities.term_vel_N_rai * If(Y.N_rai))
    @. dY.ρq_rai += -∂(aux.precip_velocities.term_vel_rai * If(Y.ρq_rai))

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.N_rai += fcc(aux.precip_velocities.term_vel_N_rai, Y.N_rai)
    @. dY.ρq_rai += fcc(aux.precip_velocities.term_vel_rai, Y.ρq_rai)

    return dY
end

"""
   Additional source terms
"""
@inline function sources_tendency!(sp::K1D.AbstractPrecipitationStyle, dY, Y, aux, t)
    error("sources_tendency not implemented for a given $sp")
end
@inline function sources_tendency!(::Union{K1D.NoPrecipitation, K1D.Precipitation0M}, dY, Y, aux, t) end
@inline function sources_tendency!(::K1D.Precipitation1M, dY, Y, aux, t)

    @. dY.ρq_tot += aux.moisture_variables.ρ * aux.precip_sources.S_q_tot
    @. dY.ρq_liq += aux.moisture_variables.ρ * aux.precip_sources.S_q_liq
    @. dY.ρq_rai += aux.moisture_variables.ρ * aux.precip_sources.S_q_rai

    return dY
end
@inline function sources_tendency!(::K1D.Precipitation2M, dY, Y, aux, t)

    @. dY.ρq_tot += aux.moisture_variables.ρ * aux.precip_sources.S_q_tot
    @. dY.ρq_liq += aux.moisture_variables.ρ * aux.precip_sources.S_q_liq
    @. dY.ρq_rai += aux.moisture_variables.ρ * aux.precip_sources.S_q_rai
    @. dY.N_liq += aux.precip_sources.S_N_liq
    @. dY.N_rai += aux.precip_sources.S_N_rai

    return dY
end
