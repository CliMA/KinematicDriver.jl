"""
   Different components of the ODE `rhs!` function,
   depending on the precipitation types.
"""

"""
     Helper functions for broadcasting and populating the auxiliary state
"""
@inline function precip_helper_sources!(ps::Precipitation1M, box_params, q_liq, q_rai, ρ, dt)

    rf = ps.rain_formation

    FT = eltype(q_liq)
    S_q_liq::FT = FT(0)
    S_q_rai::FT = FT(0)

    tmp::FT = FT(0)
    # autoconversion liquid to rain
    if rf isa CMP.Acnv1M{FT}
        tmp = CM1.conv_q_liq_to_q_rai(rf, q_liq, smooth_transition = true)
    elseif typeof(rf) in [CMP.LD2004{FT}, CMP.VarTimescaleAcnv{FT}]
        tmp = CM2.conv_q_liq_to_q_rai(rf, q_liq, ρ, ; conv_args(rf, box_params)...)
    elseif typeof(rf.acnv) in [CMP.AcnvKK2000{FT}, CMP.AcnvB1994{FT}, CMP.AcnvTC1980{FT}]
        tmp = CM2.conv_q_liq_to_q_rai(rf, q_liq, ρ, ; conv_args(rf, box_params)...)
    else
        error("Unrecognized rain formation scheme")
    end
    S_qt_rain = -min(max(FT(0), q_liq / dt), tmp)
    S_q_liq += S_qt_rain
    S_q_rai -= S_qt_rain

    # accretion cloud water + rain
    if typeof(rf) in [CMP.Acnv1M{FT}, CMP.LD2004{FT}, CMP.VarTimescaleAcnv{FT}]
        tmp = CM1.accretion(ps.liquid, ps.rain, ps.sedimentation.rain, ps.ce, q_liq, q_rai, ρ)
    elseif typeof(rf.accr) in [CMP.AccrKK2000{FT}, CMP.AccrB1994{FT}, CMP.AccrTC1980{FT}]
        tmp = CM2.accretion(accr_args(rf, q_liq, q_rai, ρ)...)
    else
        error("Unrecognized rain formation scheme")
    end
    S_qr = min(max(FT(0), q_liq / dt), tmp)
    S_q_liq -= S_qr
    S_q_rai += S_qr

    return (; S_q_liq, S_q_rai)
end
@inline function precip_helper_sources!(ps::Precipitation2M, box_params, q_liq, q_rai, N_liq, N_rai, ρ, dt)

    sb2006 = ps.rain_formation

    FT = eltype(q_liq)

    S_q_liq::FT = FT(0)
    S_q_rai::FT = FT(0)
    S_N_liq::FT = FT(0)
    S_N_rai::FT = FT(0)

    if Bool(box_params.precip_sources)
        # autoconversion liquid to rain
        tmp = CM2.autoconversion(sb2006.acnv, q_liq, q_rai, ρ, N_liq)
        S_qr = min(max(FT(0), q_liq / dt), tmp.dq_rai_dt)
        S_q_liq -= S_qr
        S_q_rai += S_qr
        S_Nr = min(max(FT(0), N_liq / 2 / dt), tmp.dN_rai_dt)
        S_N_liq -= 2 * S_Nr
        S_N_rai += S_Nr

        # accretion cloud water + rain
        tmp = CM2.accretion(sb2006, q_liq, q_rai, ρ, N_liq)
        S_qr = min(max(FT(0), q_liq / dt), tmp.dq_rai_dt)
        S_q_rai += S_qr
        S_q_liq -= S_qr
        S_N_liq += -min(max(FT(0), N_liq / dt), -tmp.dN_liq_dt)
        S_N_rai += FT(0)
    end

    if Bool(box_params.precip_sinks)
        # self_collection
        tmp_l = CM2.liquid_self_collection(sb2006.acnv, q_liq, ρ, -2 * S_Nr)
        tmp_r = CM2.rain_self_collection(sb2006.pdf, sb2006.self, q_rai, ρ, N_rai)
        S_N_liq += -min(max(FT(0), N_liq / dt), -tmp_l)
        S_N_rai += -min(max(FT(0), N_rai / dt), -tmp_r)

        # rain breakup
        tmp = CM2.rain_breakup(sb2006.pdf, sb2006.brek, q_rai, ρ, N_rai, tmp_r)
        S_N_rai += min(max(FT(0), N_rai / dt), tmp_r)
    end

    return (; S_q_liq, S_q_rai, S_N_liq, S_N_rai)
end

"""
     Precompute the auxiliary values
"""
@inline function precompute_aux_precip!(sp::AbstractPrecipitationStyle, dY, Y, aux, t)
    error("precompute_aux not implemented for a given $sp")
end
@inline function precompute_aux_precip!(ps::Precipitation1M, dY, Y, aux, t)

    ρ = BP.ρ_air(aux.box_params) / (1 - Y[1])
    tmp = precip_helper_sources!(ps, aux.box_params, Y[1], Y[2], ρ, aux.TS.dt)
    dY[1] = tmp.S_q_liq
    dY[2] = tmp.S_q_rai
end
@inline function precompute_aux_precip!(ps::Precipitation2M, dY, Y, aux, t)

    ρ = BP.ρ_air(aux.box_params) / (1 - Y[1])
    tmp = precip_helper_sources!(ps, aux.box_params, Y[1], Y[2], Y[3], Y[4], ρ, aux.TS.dt)
    dY[1] = tmp.S_q_liq
    dY[2] = tmp.S_q_rai
    dY[3] = tmp.S_N_liq
    dY[4] = tmp.S_N_rai
end
