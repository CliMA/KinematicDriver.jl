"""
   Different components of the ODE `rhs!` function,
   depending on the moisture and precipitation types.
"""

"""
    Zero out previous timestep tendencies and aux.S source terms
"""
@inline function zero_tendencies!(sm::AbstractMoistureStyle, dY, Y, aux, t)
    error("zero_tendecies not implemented for a given $sm")
end
@inline function zero_tendencies!(sp::AbstractPrecipitationStyle, dY, Y, aux, t)
    error("zero_tendecies not implemented for a given $sp")
end
@inline function zero_tendencies!(::EquilibriumMoisture, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    @. dY.ρq_tot = FT(0)
end
@inline function zero_tendencies!(::NonEquilibriumMoisture, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    @. dY.ρq_tot = FT(0)
    @. dY.ρq_liq = FT(0)
    @. dY.ρq_ice = FT(0)
end
@inline function zero_tendencies!(::Union{NoPrecipitation, Precipitation0M}, dY, Y, aux, t) end
@inline function zero_tendencies!(::Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    @. dY.ρq_rai = FT(0)
    @. dY.ρq_sno = FT(0)
end
@inline function zero_tendencies!(::Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    @. dY.ρq_rai = FT(0)
    @. dY.ρq_sno = FT(0)
    @. dY.N_liq = FT(0)
    @. dY.N_rai = FT(0)
    @. dY.N_aer = FT(0)
end

"""
     Helper functions for broadcasting and populating the auxiliary state
"""
@inline function moisture_helper_vars_eq_ρθq(thermo_params, ρq_tot, ρ, θ_liq_ice)

    ts = TD.PhaseEquil_ρθq(thermo_params, ρ, θ_liq_ice, ρq_tot / ρ)

    q_tot = TD.total_specific_humidity(thermo_params, ts)
    q_liq = TD.liquid_specific_humidity(thermo_params, ts)
    q_ice = TD.ice_specific_humidity(thermo_params, ts)

    ρ_dry = ρ - ρq_tot
    p = TD.air_pressure(thermo_params, ts)
    T = TD.air_temperature(thermo_params, ts)
    θ_dry = TD.dry_pottemp(thermo_params, T, ρ_dry)

    return (; ts, q_tot, q_liq, q_ice, ρ_dry, p, T, θ_dry)
end
@inline function moisture_helper_vars_eq_ρdTq(thermo_params, ρq_tot, ρ_dry, T)

    ρ = ρ_dry .+ ρq_tot
    ts = TD.PhaseEquil_ρTq(thermo_params, ρ, T, ρq_tot / ρ)
    p = TD.air_pressure(thermo_params, ts)

    q_tot = TD.total_specific_humidity(thermo_params, ts)
    q_liq = TD.liquid_specific_humidity(thermo_params, ts)
    q_ice = TD.ice_specific_humidity(thermo_params, ts)

    θ_liq_ice = TD.liquid_ice_pottemp(thermo_params, ts)
    θ_dry = TD.dry_pottemp(thermo_params, T, ρ_dry)

    return (; ts, q_tot, q_liq, q_ice, ρ, p, θ_liq_ice, θ_dry)
end
@inline function moisture_helper_vars_neq_ρθq(thermo_params, ρq_tot, ρq_liq, ρq_ice, ρ, θ_liq_ice)

    q_tot = ρq_tot / ρ
    q_liq = ρq_liq / ρ
    q_ice = ρq_ice / ρ
    q = TD.PhasePartition(q_tot, q_liq, q_ice)

    ts = TD.PhaseNonEquil_ρθq(thermo_params, ρ, θ_liq_ice, q)

    ρ_dry = ρ - ρq_tot
    p = TD.air_pressure(thermo_params, ts)
    T = TD.air_temperature(thermo_params, ts)
    θ_dry = TD.dry_pottemp(thermo_params, T, ρ_dry)

    return (; ts, q_tot, q_liq, q_ice, ρ_dry, p, T, θ_dry)
end
@inline function moisture_helper_vars_neq_ρdTq(thermo_params, ρq_tot, ρq_liq, ρq_ice, ρ_dry, T)

    ρ = ρ_dry .+ ρq_tot

    q_tot = ρq_tot / ρ
    q_liq = ρq_liq / ρ
    q_ice = ρq_ice / ρ
    q = TD.PhasePartition(q_tot, q_liq, q_ice)

    ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
    p = TD.air_pressure(thermo_params, ts)

    θ_liq_ice = TD.liquid_ice_pottemp(thermo_params, ts)
    θ_dry = TD.dry_pottemp(thermo_params, T, ρ_dry)

    return (; ts, q_tot, q_liq, q_ice, ρ, p, θ_liq_ice, θ_dry)
end
@inline function moisture_helper_sources(thermo_params, ne, ρ, T, q_tot, q_liq, q_ice)

    q = TD.PhasePartition(q_tot, q_liq, q_ice)

    # from sat adjst
    #ts = TD.PhaseEquil_ρθq(params, ρ, θ_liq_ice, ρq_tot / ρ)
    ts_eq = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)
    q_eq = TD.PhasePartition(thermo_params, ts_eq)

    S_q_liq = CMNe.conv_q_vap_to_q_liq_ice(ne.liquid, q_eq, q)
    S_q_ice = CMNe.conv_q_vap_to_q_liq_ice(ne.ice, q_eq, q)

    return (; S_q_liq, S_q_ice)
end
@inline function precip_helper_vars(ρq_rai, ρq_sno, ρ)

    q_rai = ρq_rai / ρ
    q_sno = ρq_sno / ρ

    return (; q_rai, q_sno)
end
@inline function precip_helper_sources!(ps::Precipitation0M, thermo_params, ts, q_tot, q_liq, q_ice, dt)

    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    qsat = TD.q_vap_saturation(thermo_params, ts)
    λ = TD.liquid_fraction(thermo_params, ts)

    S_qt = -min(max(0, (q.liq + q.ice) / dt), -CM0.remove_precipitation(ps.params, q, qsat))

    S_q_rai = -S_qt * λ
    S_q_sno = -S_qt * (1 - λ)
    S_q_tot = S_qt
    S_q_liq = S_qt * λ
    S_q_ice = S_qt * (1 - λ)
    #θ_liq_ice_tendency -= S_qt / Π_m / c_pm * (L_v0 * λ + L_s0 * (1 - λ))

    return (; S_q_tot, S_q_liq, S_q_ice, S_q_rai, S_q_sno)
end
@inline function precip_helper_sources!(
    ps::Precipitation1M,
    common_params,
    thermo_params,
    air_params,
    ts,
    q_tot,
    q_liq,
    q_ice,
    q_rai,
    q_sno,
    T,
    ρ,
    dt,
)

    FT = eltype(q_tot)

    S_q_rai::FT = FT(0)
    S_q_sno::FT = FT(0)
    S_q_tot::FT = FT(0)
    S_q_liq::FT = FT(0)
    S_q_ice::FT = FT(0)
    S_q_vap::FT = FT(0) # added for testing

    T_fr::FT = TD.Parameters.T_freeze(thermo_params)
    c_vl::FT = TD.Parameters.cv_l(thermo_params)
    c_vm::FT = TD.cv_m(thermo_params, ts)
    Rm::FT = TD.gas_constant_air(thermo_params, ts)
    Lf::FT = TD.latent_heat_fusion(thermo_params, ts)

    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)

    if Bool(common_params.precip_sources)
        tmp::FT = FT(0)
        rf = ps.rain_formation
        # autoconversion liquid to rain and ice to snow
        # TODO - can we do it in a more elegant way?
        if rf isa CMP.Acnv1M{FT}
            tmp = CM1.conv_q_liq_to_q_rai(rf, q.liq, smooth_transition = true)
        elseif typeof(rf) in [CMP.LD2004{FT}, CMP.VarTimescaleAcnv{FT}]
            tmp = CM2.conv_q_liq_to_q_rai(rf, q.liq, ρ; N_d = common_params.prescribed_Nd)
        elseif typeof(rf.acnv) in [CMP.AcnvKK2000{FT}, CMP.AcnvB1994{FT}, CMP.AcnvTC1980{FT}]
            tmp = CM2.conv_q_liq_to_q_rai(rf, q.liq, ρ; N_d = common_params.prescribed_Nd)
        else
            error("Unrecognized rain formation scheme")
        end
        S_qt_rain = -min(max(FT(0), q.liq / dt), tmp)
        S_qt_snow = -min(max(FT(0), q.ice / dt), CM1.conv_q_ice_to_q_sno(ps.ice, air_params, thermo_params, q, ρ, T))
        S_q_rai -= S_qt_rain
        S_q_sno -= S_qt_snow
        S_q_tot += S_qt_rain + S_qt_snow
        S_q_liq += S_qt_rain
        S_q_ice += S_qt_snow
        #θ_liq_ice_tendency -= 1 / Π_m / c_pm * (L_v0 * S_qt_rain + L_s0 * S_qt_snow)

        # accretion cloud water + rain
        if typeof(rf) in [CMP.Acnv1M{FT}, CMP.LD2004{FT}, CMP.VarTimescaleAcnv{FT}]
            tmp = CM1.accretion(ps.liquid, ps.rain, ps.sedimentation.rain, ps.ce, q.liq, q_rai, ρ)
        elseif typeof(rf.accr) in [CMP.AccrKK2000{FT}, CMP.AccrB1994{FT}]
            tmp = CM2.accretion(rf, q.liq, q_rai, ρ)
        elseif rf.accr isa CMP.AccrTC1980{FT}
            tmp = CM2.accretion(rf, q.liq, q_rai)
        else
            error("Unrecognized rain formation scheme")
        end
        S_qr = min(max(FT(0), q.liq / dt), tmp)
        S_q_rai += S_qr
        S_q_tot -= S_qr
        S_q_liq -= S_qr
        #θ_liq_ice_tendency += S_qr / Π_m / c_pm * L_v0

        # accretion cloud ice + snow
        S_qs =
            min(max(FT(0), q.ice / dt), CM1.accretion(ps.ice, ps.snow, ps.sedimentation.snow, ps.ce, q.ice, q_sno, ρ))
        S_q_sno += S_qs
        S_q_tot -= S_qs
        S_q_ice -= S_qs
        #θ_liq_ice_tendency += S_qs / Π_m / c_pm * L_s0

        # sink of cloud water via accretion cloud water + snow
        S_qt =
            -min(
                max(FT(0), q.liq / dt),
                CM1.accretion(ps.liquid, ps.snow, ps.sedimentation.snow, ps.ce, q.liq, q_sno, ρ),
            )
        if T < T_fr # cloud droplets freeze to become snow)
            S_q_sno -= S_qt
            S_q_tot += S_qt
            S_q_liq += S_qt
            #θ_liq_ice_tendency -= S_qt / Π_m / c_pm * Lf * (1 + Rm / c_vm)
        else # snow melts, both cloud water and snow become rain
            α::FT = c_vl / Lf * (T - T_fr)
            S_q_tot += S_qt
            S_q_liq += S_qt
            S_q_sno += S_qt * α
            S_q_rai -= S_qt * (1 + α)
            #θ_liq_ice_tendency += S_qt / Π_m / c_pm * (Lf * (1 + Rm / c_vm) * α - L_v0)
        end

        # sink of cloud ice via accretion cloud ice - rain
        S_qt =
            -min(max(FT(0), q.ice / dt), CM1.accretion(ps.ice, ps.rain, ps.sedimentation.rain, ps.ce, q.ice, q_rai, ρ))
        # sink of rain via accretion cloud ice - rain
        S_qr =
            -min(
                max(FT(0), q_rai / dt),
                CM1.accretion_rain_sink(ps.rain, ps.ice, ps.sedimentation.rain, ps.ce, q.ice, q_rai, ρ),
            )
        S_q_tot += S_qt
        S_q_ice += S_qt
        S_q_rai += S_qr
        S_q_sno += -(S_qt + S_qr)
        #θ_liq_ice_tendency -= 1 / Π_m / c_pm * (S_qr * Lf * (1 + Rm / c_vm) + S_qt * L_s0)

        # accretion rain - snow
        if T < T_fr
            S_qs = min(
                max(FT(0), q_rai / dt),
                CM1.accretion_snow_rain(
                    ps.snow,
                    ps.rain,
                    ps.sedimentation.snow,
                    ps.sedimentation.rain,
                    ps.ce,
                    q_sno,
                    q_rai,
                    ρ,
                ),
            )
        else
            S_qs =
                -min(
                    max(FT(0), q_sno / dt),
                    CM1.accretion_snow_rain(
                        ps.rain,
                        ps.snow,
                        ps.sedimentation.rain,
                        ps.sedimentation.snow,
                        ps.ce,
                        q_rai,
                        q_sno,
                        ρ,
                    ),
                )
        end
        S_q_sno += S_qs
        S_q_rai -= S_qs
        #θ_liq_ice_tendency += S_qs * Lf / Π_m / c_vm
    end

    if Bool(common_params.precip_sinks)
        # evaporation
        S_qr =
            -min(
                max(FT(0), q_rai / dt),
                -CM1.evaporation_sublimation(ps.rain, ps.sedimentation.rain, air_params, thermo_params, q, q_rai, ρ, T),
            )
        S_q_rai += S_qr
        S_q_tot -= S_qr
        S_q_vap -= S_qr

        # melting
        S_qs =
            -min(
                max(FT(0), q_sno / dt),
                CM1.snow_melt(ps.snow, ps.sedimentation.snow, air_params, thermo_params, q_sno, ρ, T),
            )
        S_q_rai -= S_qs
        S_q_sno += S_qs

        # deposition, sublimation
        tmp = CM1.evaporation_sublimation(ps.snow, ps.sedimentation.snow, air_params, thermo_params, q, q_sno, ρ, T)
        if tmp > 0
            S_qs = min(max(FT(0), q_vap / dt), tmp)
        else
            S_qs = -min(max(FT(0), q_sno / dt), -tmp)
        end
        S_q_sno += S_qs
        S_q_tot -= S_qs
        S_q_vap -= S_qs

        #θ_liq_ice_tendency +=
        #    1 / Π_m / c_pm * (
        #        S_qr_evap * (L_v - R_v * T) * (1 + R_m / c_vm) +
        #        S_qs_sub_dep * (L_s - R_v * T) * (1 + R_m / c_vm) +
        #        S_qs_melt * L_f * (1 + R_m / c_vm)
        #)
    end

    term_vel_rai = CM1.terminal_velocity(ps.rain, ps.sedimentation.rain, ρ, q_rai)
    term_vel_sno = CM1.terminal_velocity(ps.snow, ps.sedimentation.snow, ρ, q_sno)

    return (; S_q_tot, S_q_liq, S_q_ice, S_q_rai, S_q_sno, S_q_vap, term_vel_rai, term_vel_sno)
end
@inline function precip_helper_sources!(
    ps::Precipitation2M,
    common_params,
    thermo_params,
    air_params,
    q_tot,
    q_liq,
    q_rai,
    N_liq,
    N_rai,
    T,
    ρ,
    dt,
)

    sb2006 = ps.rain_formation
    vel_scheme = ps.sedimentation

    FT = eltype(q_tot)

    S_q_rai::FT = FT(0)
    S_q_tot::FT = FT(0)
    S_q_liq::FT = FT(0)
    S_q_vap::FT = FT(0) # added for testing
    S_N_liq::FT = FT(0)
    S_N_rai::FT = FT(0)

    q = TD.PhasePartition(q_tot, q_liq, FT(0))

    if Bool(common_params.precip_sources)
        # autoconversion liquid to rain
        tmp = CM2.autoconversion(sb2006.acnv, q.liq, q_rai, ρ, N_liq)
        S_qr = min(max(FT(0), q.liq / dt), tmp.dq_rai_dt)
        S_q_rai += S_qr
        S_q_tot -= S_qr
        S_q_liq -= S_qr
        S_Nr = min(max(FT(0), N_liq / 2 / dt), tmp.dN_rai_dt)
        S_N_rai += S_Nr
        S_N_liq -= 2 * S_Nr

        # self_collection
        tmp_l = CM2.liquid_self_collection(sb2006.acnv, q.liq, ρ, -2 * S_Nr)
        tmp_r = CM2.rain_self_collection(sb2006.pdf, sb2006.self, q_rai, ρ, N_rai)
        S_N_liq += -min(max(FT(0), N_liq / dt), -tmp_l)
        S_N_rai += -min(max(FT(0), N_rai / dt), -tmp_r)

        # rain breakup
        tmp = CM2.rain_breakup(sb2006.pdf, sb2006.brek, q_rai, ρ, N_rai, tmp_r)
        S_N_rai += min(max(FT(0), N_rai / dt), tmp_r)

        # accretion cloud water + rain
        tmp = CM2.accretion(sb2006, q.liq, q_rai, ρ, N_liq)
        S_qr = min(max(FT(0), q.liq / dt), tmp.dq_rai_dt)
        S_q_rai += S_qr
        S_q_tot -= S_qr
        S_q_liq -= S_qr
        S_N_liq += -min(max(FT(0), N_liq / dt), -tmp.dN_liq_dt)
        S_N_rai += FT(0)
    end

    if Bool(common_params.precip_sinks)
        # evaporation
        tmp = CM2.rain_evaporation(sb2006, air_params, thermo_params, q, q_rai, ρ, N_rai, T)
        S_Nr = -min(max(FT(0), N_rai / dt), -tmp[1])
        S_qr = -min(max(FT(0), q_rai / dt), -tmp[2])
        S_q_rai += S_qr
        S_q_tot -= S_qr
        S_q_vap -= S_qr
        S_N_rai += S_Nr
    end

    vt = CM2.rain_terminal_velocity(sb2006, vel_scheme, q_rai, ρ, N_rai)
    term_vel_N_rai = vt[1]
    term_vel_rai = vt[2]

    return (; S_q_tot, S_q_liq, S_q_rai, S_N_liq, S_N_rai, S_q_vap, term_vel_rai, term_vel_N_rai)
end

"""
     Precompute the auxiliary values
"""
@inline function precompute_aux_thermo!(sm::AbstractMoistureStyle, dY, Y, aux, t)
    error("precompute_aux not implemented for a given $sm")
end
@inline function precompute_aux_thermo!(::EquilibriumMoisture_ρθq, dY, Y, aux, t)
    
    tmp = @. moisture_helper_vars_eq_ρθq(
        aux.thermo_params,
        Y.ρq_tot,
        aux.moisture_variables.ρ,
        aux.moisture_variables.θ_liq_ice,
    )
    aux.moisture_variables.ρ_dry = tmp.ρ_dry
    aux.moisture_variables.p = tmp.p
    aux.moisture_variables.T = tmp.T
    aux.moisture_variables.θ_dry = tmp.θ_dry
    aux.moisture_variables.q_tot = tmp.q_tot
    aux.moisture_variables.q_liq = tmp.q_liq
    aux.moisture_variables.q_ice = tmp.q_ice
    aux.moisture_variables.ts = tmp.ts

end
@inline function precompute_aux_thermo!(::EquilibriumMoisture_ρdTq, dY, Y, aux, t)

    tmp = @. moisture_helper_vars_eq_ρdTq(
        aux.thermo_params,
        Y.ρq_tot,
        aux.moisture_variables.ρ_dry,
        aux.moisture_variables.T,
    )
    aux.moisture_variables.ρ = tmp.ρ
    aux.moisture_variables.p = tmp.p
    aux.moisture_variables.θ_liq_ice = tmp.θ_liq_ice
    aux.moisture_variables.θ_dry = tmp.θ_dry
    aux.moisture_variables.q_tot = tmp.q_tot
    aux.moisture_variables.q_liq = tmp.q_liq
    aux.moisture_variables.q_ice = tmp.q_ice
    aux.moisture_variables.ts = tmp.ts

end
@inline function precompute_aux_thermo!(::NonEquilibriumMoisture_ρθq, dY, Y, aux, t)
    
    tmp = @. moisture_helper_vars_neq_ρθq(
        aux.thermo_params,
        Y.ρq_tot,
        Y.ρq_liq,
        Y.ρq_ice,
        aux.moisture_variables.ρ,
        aux.moisture_variables.θ_liq_ice,
    )
    aux.moisture_variables.ρ_dry = tmp.ρ_dry
    aux.moisture_variables.p = tmp.p
    aux.moisture_variables.T = tmp.T
    aux.moisture_variables.θ_dry = tmp.θ_dry
    aux.moisture_variables.q_tot = tmp.q_tot
    aux.moisture_variables.q_liq = tmp.q_liq
    aux.moisture_variables.q_ice = tmp.q_ice
    aux.moisture_variables.ts = tmp.ts

end
@inline function precompute_aux_thermo!(::NonEquilibriumMoisture_ρdTq, dY, Y, aux, t)
    
    tmp = @. moisture_helper_vars_neq_ρdTq(
        aux.thermo_params,
        Y.ρq_tot,
        Y.ρq_liq,
        Y.ρq_ice,
        aux.moisture_variables.ρ_dry,
        aux.moisture_variables.T,
    )
    aux.moisture_variables.ρ = tmp.ρ
    aux.moisture_variables.p = tmp.p
    aux.moisture_variables.θ_liq_ice = tmp.θ_liq_ice
    aux.moisture_variables.θ_dry = tmp.θ_dry
    aux.moisture_variables.q_tot = tmp.q_tot
    aux.moisture_variables.q_liq = tmp.q_liq
    aux.moisture_variables.q_ice = tmp.q_ice
    aux.moisture_variables.ts = tmp.ts

end
@inline function precompute_aux_moisture_sources!(sm::AbstractMoistureStyle, dY, Y, aux, t)
    error("precompute_aux not implemented for a given $sm")
end
@inline function precompute_aux_moisture_sources!(::EquilibriumMoisture, dY, Y, aux, t) end
@inline function precompute_aux_moisture_sources!(ne::NonEquilibriumMoisture, dY, Y, aux, t)
    tmp = @. moisture_helper_sources(
        aux.thermo_params,
        ne,
        aux.moisture_variables.ρ,
        aux.moisture_variables.T,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.moisture_variables.q_ice,
    )
    aux.moisture_sources.S_q_liq = tmp.S_q_liq
    aux.moisture_sources.S_q_ice = tmp.S_q_ice

end
@inline function precompute_aux_precip!(sp::AbstractPrecipitationStyle, dY, Y, aux, t)
    error("precompute_aux not implemented for a given $sp")
end
@inline function precompute_aux_precip!(::NoPrecipitation, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    aux.precip_variables.q_rai = FT(0)
    aux.precip_variables.q_sno = FT(0)

    aux.precip_sources.S_q_rai = FT(0)
    aux.precip_sources.S_q_sno = FT(0)
    aux.precip_sources.S_q_tot = FT(0)
    aux.precip_sources.S_q_liq = FT(0)
    aux.precip_sources.S_q_ice = FT(0)
end
@inline function precompute_aux_precip!(ps::Precipitation0M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    aux.precip_variables.q_rai = FT(0)
    aux.precip_variables.q_sno = FT(0)

    tmp = @. precip_helper_sources!(
        ps,
        aux.thermo_params,
        aux.moisture_variables.ts,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.moisture_variables.q_ice,
        aux.TS.dt,
    )
    aux.precip_sources.S_q_rai = tmp.S_q_rai
    aux.precip_sources.S_q_sno = tmp.S_q_sno
    aux.precip_sources.S_q_tot = tmp.S_q_tot
    aux.precip_sources.S_q_liq = tmp.S_q_liq
    aux.precip_sources.S_q_ice = tmp.S_q_ice
end
@inline function precompute_aux_precip!(ps::Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_rai)

    tmp_v = @. precip_helper_vars(Y.ρq_rai, Y.ρq_sno, aux.moisture_variables.ρ)
    aux.precip_variables.q_rai = tmp_v.q_rai
    aux.precip_variables.q_sno = tmp_v.q_sno

    tmp = @. precip_helper_sources!(
        ps,
        aux.common_params,
        aux.thermo_params,
        aux.air_params,
        aux.moisture_variables.ts,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.moisture_variables.q_ice,
        aux.precip_variables.q_rai,
        aux.precip_variables.q_sno,
        aux.moisture_variables.T,
        aux.moisture_variables.ρ,
        aux.TS.dt,
    )
    aux.precip_sources.S_q_rai = tmp.S_q_rai
    aux.precip_sources.S_q_sno = tmp.S_q_sno
    aux.precip_sources.S_q_tot = tmp.S_q_tot
    aux.precip_sources.S_q_liq = tmp.S_q_liq
    aux.precip_sources.S_q_ice = tmp.S_q_ice

    aux.precip_velocities.term_vel_rai = tmp.term_vel_rai
    aux.precip_velocities.term_vel_sno = tmp.term_vel_sno
end
@inline function precompute_aux_precip!(ps::Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_rai)

    aux.precip_variables.q_rai = Y.ρq_rai ./ aux.moisture_variables.ρ
    aux.precip_variables.N_liq = Y.N_liq
    aux.precip_variables.N_rai = Y.N_rai
    tmp = @. precip_helper_sources!(
        ps,
        aux.common_params,
        aux.thermo_params,
        aux.air_params,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.precip_variables.q_rai,
        aux.precip_variables.N_liq,
        aux.precip_variables.N_rai,
        aux.moisture_variables.T,
        aux.moisture_variables.ρ,
        aux.TS.dt,
    )
    aux.precip_sources.S_q_rai = tmp.S_q_rai
    aux.precip_sources.S_q_tot = tmp.S_q_tot
    aux.precip_sources.S_q_liq = tmp.S_q_liq
    aux.precip_sources.S_N_liq = tmp.S_N_liq
    aux.precip_sources.S_N_rai = tmp.S_N_rai

    aux.precip_velocities.term_vel_N_rai = tmp.term_vel_N_rai
    aux.precip_velocities.term_vel_rai = tmp.term_vel_rai
end

"""
   Additional source terms
"""
@inline function sources_tendency!(sm::AbstractMoistureStyle, dY, Y, aux, t)
    error("sources_tendency not implemented for a given $sm")
end
@inline function sources_tendency!(sp::AbstractPrecipitationStyle, dY, Y, aux, t)
    error("sources_tendency not implemented for a given $sp")
end
@inline function sources_tendency!(::EquilibriumMoisture, dY, Y, aux, t)

    @. dY.ρq_tot += aux.moisture_variables.ρ * aux.precip_sources.S_q_tot

    return dY
end
@inline function sources_tendency!(::NonEquilibriumMoisture, dY, Y, aux, t)

    @. dY.ρq_liq += aux.moisture_variables.ρ * aux.moisture_sources.S_q_liq
    @. dY.ρq_ice += aux.moisture_variables.ρ * aux.moisture_sources.S_q_ice

    @. dY.ρq_tot += aux.moisture_variables.ρ * aux.precip_sources.S_q_tot
    @. dY.ρq_liq += aux.moisture_variables.ρ * aux.precip_sources.S_q_liq
    @. dY.ρq_ice += aux.moisture_variables.ρ * aux.precip_sources.S_q_ice

    return dY
end
@inline function sources_tendency!(::Union{NoPrecipitation, Precipitation0M}, dY, Y, aux, t) end
@inline function sources_tendency!(::Precipitation1M, dY, Y, aux, t)

    @. dY.ρq_rai += aux.moisture_variables.ρ * aux.precip_sources.S_q_rai
    @. dY.ρq_sno += aux.moisture_variables.ρ * aux.precip_sources.S_q_sno

    return dY
end
@inline function sources_tendency!(::Precipitation2M, dY, Y, aux, t)

    @. dY.ρq_rai += aux.moisture_variables.ρ * aux.precip_sources.S_q_rai
    @. dY.ρq_sno += aux.moisture_variables.ρ * aux.precip_sources.S_q_sno
    @. dY.N_liq += aux.precip_sources.S_N_liq
    @. dY.N_rai += aux.precip_sources.S_N_rai

    @. dY.N_liq += aux.activation_sources.S_N_liq
    @. dY.N_aer += aux.activation_sources.S_N_aer

    return dY
end
