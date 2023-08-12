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
    @. dY.N_liq = FT(0)
    @. dY.N_rai = FT(0)
    @. dY.N_aer = FT(0)
end

"""
     Helper functions for broadcasting and populating the auxiliary state
"""
@inline function moisture_helper_vars_eq_ρθq(params, ρq_tot, ρ, θ_liq_ice)

    thermo_params = KP.thermodynamics_params(params)

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
@inline function moisture_helper_vars_eq_ρdTq(params, ρq_tot, ρ_dry, T)

    thermo_params = KP.thermodynamics_params(params)

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
@inline function moisture_helper_vars_neq_ρθq(params, ρq_tot, ρq_liq, ρq_ice, ρ, θ_liq_ice)

    thermo_params = KP.thermodynamics_params(params)

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
@inline function moisture_helper_vars_neq_ρdTq(params, ρq_tot, ρq_liq, ρq_ice, ρ_dry, T)

    thermo_params = KP.thermodynamics_params(params)

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
@inline function moisture_helper_sources(params, ρ, T, q_tot, q_liq, q_ice)

    thermo_params = KP.thermodynamics_params(params)
    microphys_params = KP.microphysics_params(params)

    q = TD.PhasePartition(q_tot, q_liq, q_ice)

    # from sat adjst
    #ts = TD.PhaseEquil_ρθq(params, ρ, θ_liq_ice, ρq_tot / ρ)
    ts_eq = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)
    q_eq = TD.PhasePartition(thermo_params, ts_eq)

    S_q_liq = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, CMT.LiquidType(), q_eq, q)
    S_q_ice = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, CMT.IceType(), q_eq, q)

    return (; S_q_liq, S_q_ice)
end
@inline function precip_helper_vars(ρq_rai, ρq_sno, ρ)

    q_rai = ρq_rai / ρ
    q_sno = ρq_sno / ρ

    return (; q_rai, q_sno)
end
@inline function precip_helper_sources_0M!(params, ts, q_tot, q_liq, q_ice, dt)

    thermo_params = KP.thermodynamics_params(params)
    microphys_params = KP.microphysics_params(params)

    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    qsat = TD.q_vap_saturation(thermo_params, ts)
    λ = TD.liquid_fraction(thermo_params, ts)

    S_qt = -min(max(0, (q.liq + q.ice) / dt), -CM0.remove_precipitation(microphys_params, q, qsat))

    S_q_rai = -S_qt * λ
    S_q_sno = -S_qt * (1 - λ)
    S_q_tot = S_qt
    S_q_liq = S_qt * λ
    S_q_ice = S_qt * (1 - λ)
    #θ_liq_ice_tendency -= S_qt / Π_m / c_pm * (L_v0 * λ + L_s0 * (1 - λ))

    return (; S_q_tot, S_q_liq, S_q_ice, S_q_rai, S_q_sno)
end
@inline function precip_helper_sources_1M!(params, ts, q_tot, q_liq, q_ice, q_rai, q_sno, T, ρ, dt, ps)

    microphys_params = KP.microphysics_params(params)
    thermo_params = KP.thermodynamics_params(params)

    rf = ps.rain_formation

    FT = eltype(q_tot)

    S_q_rai::FT = FT(0)
    S_q_sno::FT = FT(0)
    S_q_tot::FT = FT(0)
    S_q_liq::FT = FT(0)
    S_q_ice::FT = FT(0)
    S_q_vap::FT = FT(0) # added for testing

    T_fr::FT = KP.T_freeze(params)
    c_vl::FT = KP.cv_l(params)
    c_vm::FT = TD.cv_m(thermo_params, ts)
    Rm::FT = TD.gas_constant_air(thermo_params, ts)
    Lf::FT = TD.latent_heat_fusion(thermo_params, ts)

    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)

    if Bool(params.precip_sources)
        tmp::FT = FT(0)
        # autoconversion liquid to rain and ice to snow
        if rf isa OneMomentRainFormation
            tmp = CM1.conv_q_liq_to_q_rai(microphys_params, q.liq, smooth_transition = true)
        elseif typeof(rf) in [CMT.KK2000Type, CMT.B1994Type, CMT.TC1980Type, CMT.LD2004Type]
            tmp = CM2.conv_q_liq_to_q_rai(microphys_params, rf, q.liq, ρ, ; conv_args(rf, params)...)
        else
            error("Unrecognized rain formation scheme")
        end
        S_qt_rain = -min(max(0, q.liq / dt), tmp)
        S_qt_snow = -min(max(0, q.ice / dt), CM1.conv_q_ice_to_q_sno(microphys_params, q, ρ, T))
        S_q_rai -= S_qt_rain
        S_q_sno -= S_qt_snow
        S_q_tot += S_qt_rain + S_qt_snow
        S_q_liq += S_qt_rain
        S_q_ice += S_qt_snow
        #θ_liq_ice_tendency -= 1 / Π_m / c_pm * (L_v0 * S_qt_rain + L_s0 * S_qt_snow)

        if rf isa OneMomentRainFormation || rf isa CMT.LD2004Type
            tmp = CM1.accretion(microphys_params, CMT.LiquidType(), CMT.RainType(), q.liq, q_rai, ρ)
        elseif typeof(rf) in [CMT.KK2000Type, CMT.B1994Type, CMT.TC1980Type]
            tmp = CM2.accretion(microphys_params, accr_args(rf, q.liq, q_rai, ρ)...)
        else
            error("Unrecognized rain formation scheme")
        end

        # accretion cloud water + rain
        S_qr = min(max(0, q.liq / dt), tmp)
        S_q_rai += S_qr
        S_q_tot -= S_qr
        S_q_liq -= S_qr
        #θ_liq_ice_tendency += S_qr / Π_m / c_pm * L_v0

        # accretion cloud ice + snow
        S_qs = min(max(0, q.ice / dt), CM1.accretion(microphys_params, CMT.IceType(), CMT.SnowType(), q.ice, q_sno, ρ))
        S_q_sno += S_qs
        S_q_tot -= S_qs
        S_q_ice -= S_qs
        #θ_liq_ice_tendency += S_qs / Π_m / c_pm * L_s0

        # sink of cloud water via accretion cloud water + snow
        S_qt =
            -min(max(0, q.liq / dt), CM1.accretion(microphys_params, CMT.LiquidType(), CMT.SnowType(), q.liq, q_sno, ρ))
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
        S_qt = -min(max(0, q.ice / dt), CM1.accretion(microphys_params, CMT.IceType(), CMT.RainType(), q.ice, q_rai, ρ))
        # sink of rain via accretion cloud ice - rain
        S_qr = -min(max(0, q_rai / dt), CM1.accretion_rain_sink(microphys_params, q.ice, q_rai, ρ))
        S_q_tot += S_qt
        S_q_ice += S_qt
        S_q_rai += S_qr
        S_q_sno += -(S_qt + S_qr)
        #θ_liq_ice_tendency -= 1 / Π_m / c_pm * (S_qr * Lf * (1 + Rm / c_vm) + S_qt * L_s0)

        # accretion rain - snow
        if T < T_fr
            S_qs = min(
                max(0, q_rai / dt),
                CM1.accretion_snow_rain(microphys_params, CMT.SnowType(), CMT.RainType(), q_sno, q_rai, ρ),
            )
        else
            S_qs =
                -min(
                    max(0, q_sno / dt),
                    CM1.accretion_snow_rain(microphys_params, CMT.RainType(), CMT.SnowType(), q_rai, q_sno, ρ),
                )
        end
        S_q_sno += S_qs
        S_q_rai -= S_qs
        #θ_liq_ice_tendency += S_qs * Lf / Π_m / c_vm
    end

    if Bool(params.precip_sinks)
        # evaporation
        S_qr = -min(max(0, q_rai / dt), -CM1.evaporation_sublimation(microphys_params, CMT.RainType(), q, q_rai, ρ, T))
        S_q_rai += S_qr
        S_q_tot -= S_qr
        S_q_vap -= S_qr

        # melting
        S_qs = -min(max(0, q_sno / dt), CM1.snow_melt(microphys_params, q_sno, ρ, T))
        S_q_rai -= S_qs
        S_q_sno += S_qs

        # deposition, sublimation
        tmp = CM1.evaporation_sublimation(microphys_params, CMT.SnowType(), q, q_sno, ρ, T)
        if tmp > 0
            S_qs = min(max(0, q_vap / dt), tmp)
        else
            S_qs = -min(max(0, q_sno / dt), -tmp)
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

    term_vel_rai = CM1.terminal_velocity(microphys_params, CMT.RainType(), true, CMT.Chen2022Type(), ρ, q_rai)
    term_vel_sno = CM1.terminal_velocity(microphys_params, CMT.SnowType(), false, CMT.Chen2022Type(), ρ, q_sno)

    # term_vel_rai = CM1.terminal_velocity(microphys_params, CMT.RainType(), ρ, q_rai)
    # term_vel_sno = CM1.terminal_velocity(microphys_params, CMT.SnowType(), ρ, q_sno)

    return (; S_q_tot, S_q_liq, S_q_ice, S_q_rai, S_q_sno, S_q_vap, term_vel_rai, term_vel_sno)
end
@inline function precip_helper_sources_2M!(params, q_tot, q_liq, q_rai, N_liq, N_rai, T, ρ, dt, ps)

    microphys_params = KP.microphysics_params(params)

    rf = ps.rain_formation

    FT = eltype(q_tot)

    S_q_rai::FT = FT(0)
    S_q_tot::FT = FT(0)
    S_q_liq::FT = FT(0)
    S_q_vap::FT = FT(0) # added for testing
    S_N_liq::FT = FT(0)
    S_N_rai::FT = FT(0)

    q = TD.PhasePartition(q_tot, q_liq, FT(0))

    # autoconversion liquid to rain
    tmp = CM2.autoconversion(microphys_params, rf, q.liq, q_rai, ρ, N_liq)
    S_qr = min(max(0, q.liq / dt), tmp.dq_rai_dt)
    S_q_rai += S_qr
    S_q_tot -= S_qr
    S_q_liq -= S_qr
    S_Nr = min(max(0, N_liq / 2 / dt), tmp.dN_rai_dt)
    S_N_rai += S_Nr
    S_N_liq -= 2 * S_Nr

    # self_collection
    tmp_l = CM2.liquid_self_collection(microphys_params, rf, q.liq, ρ, -S_qr)
    tmp_r = CM2.rain_self_collection(microphys_params, rf, q_rai, ρ, N_rai)
    S_N_liq += -min(max(0, N_liq / dt), -tmp_l)
    S_N_rai += -min(max(0, N_rai / dt), -tmp_r)

    # rain breakup
    tmp = CM2.rain_breakup(microphys_params, rf, q_rai, ρ, N_rai, tmp_r)
    S_N_rai += min(max(0, N_rai / dt), tmp_r)

    # accretion cloud water + rain
    tmp = CM2.accretion(microphys_params, rf, q.liq, q_rai, ρ, N_liq)
    S_qr = min(max(0, q.liq / dt), tmp.dq_rai_dt)
    S_q_rai += S_qr
    S_q_tot -= S_qr
    S_q_liq -= S_qr
    S_N_liq += -min(max(0, N_liq / dt), -tmp.dN_liq_dt)
    S_N_rai += FT(0)

    # evaporation
    tmp = CM2.rain_evaporation(microphys_params, rf, q, q_rai, ρ, N_rai, T)
    S_Nr = -min(max(0, N_rai / dt), -tmp[1])
    S_qr = -min(max(0, q_rai / dt), -tmp[2])
    S_q_rai += S_qr
    S_q_tot -= S_qr
    S_q_vap -= S_qr
    S_N_rai += S_Nr

    println("Chen now")
    #vt = CM2.rain_terminal_velocity(microphys_params, CMT.SB2006Type(), q_rai, ρ, N_rai)
    vt = CM2.rain_terminal_velocity(microphys_params, CMT.Chen2022Type(), q_rai, ρ, N_rai)
    term_vel_N_rai = vt[1]
    term_vel_rai = vt[2]

    return (; S_q_tot, S_q_liq, S_q_rai, S_N_liq, S_N_rai, S_q_vap, term_vel_rai, term_vel_N_rai)
end

"""
     Precompute the auxiliary values
"""
@inline function precompute_aux_prescribed_velocity!(aux, t)

    FT = eltype(aux.moisture_variables.q_tot)
    ρw = FT(ρw_helper(t, aux.params.w1, aux.params.t1))

    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw

end
@inline function precompute_aux_thermo!(sm::AbstractMoistureStyle, dY, Y, aux, t)
    error("precompute_aux not implemented for a given $sm")
end
@inline function precompute_aux_thermo!(::EquilibriumMoisture_ρθq, dY, Y, aux, t)

    # TODO - why this is not working?
    #@. aux.moisture_variables = moisture_helper_vars_eq(aux.params, Y.ρq_tot, aux.moisture_variables.ρ, aux.moisture_variables.θ_liq_ice)

    tmp =
        @. moisture_helper_vars_eq_ρθq(aux.params, Y.ρq_tot, aux.moisture_variables.ρ, aux.moisture_variables.θ_liq_ice)
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

    tmp = @. moisture_helper_vars_eq_ρdTq(aux.params, Y.ρq_tot, aux.moisture_variables.ρ_dry, aux.moisture_variables.T)
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
    #@. aux.moisture_variables = moisture_helper_vars_neq(aux.params, Y.ρq_tot, Y.ρq_liq, Y.ρq_ice, aux.moisture_variables.ρ, aux.moisture_variables.θ_liq_ice)
    #@. aux.moisture_sources = moisture_helper_sources(aux.params, aux.moisture_variables.ρ, aux.moisture_variables.T, aux.moisture_variables.q_tot, aux.moisture_variables.q_liq, aux.moisture_variables.q_ice)
    tmp = @. moisture_helper_vars_neq_ρθq(
        aux.params,
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

    tmp_s = @. moisture_helper_sources(
        aux.params,
        aux.moisture_variables.ρ,
        aux.moisture_variables.T,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.moisture_variables.q_ice,
    )
    aux.moisture_sources.S_q_liq = tmp_s.S_q_liq
    aux.moisture_sources.S_q_ice = tmp_s.S_q_ice

end
@inline function precompute_aux_thermo!(::NonEquilibriumMoisture_ρdTq, dY, Y, aux, t)
    tmp = @. moisture_helper_vars_neq_ρdTq(
        aux.params,
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

    tmp_s = @. moisture_helper_sources(
        aux.params,
        aux.moisture_variables.ρ,
        aux.moisture_variables.T,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.moisture_variables.q_ice,
    )
    aux.moisture_sources.S_q_liq = tmp_s.S_q_liq
    aux.moisture_sources.S_q_ice = tmp_s.S_q_ice

end
@inline function precompute_aux_precip!(sp::AbstractPrecipitationStyle, dY, Y, aux, t)
    error("precompute_aux not implemented for a given $sp")
end
@inline function precompute_aux_precip!(::NoPrecipitation, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    #@. aux.precip_variables = (; q_rai = FT(0), q_sno = FT(0))
    #@. aux.precip_sources = (; S_q_tot = FT(0), S_q_liq = FT(0), S_q_ice = FT(0), S_q_rai = FT(0), S_q_sno = FT(0))

    aux.precip_variables.q_rai = FT(0)
    aux.precip_variables.q_sno = FT(0)

    aux.precip_sources.S_q_rai = FT(0)
    aux.precip_sources.S_q_sno = FT(0)
    aux.precip_sources.S_q_tot = FT(0)
    aux.precip_sources.S_q_liq = FT(0)
    aux.precip_sources.S_q_ice = FT(0)
end
@inline function precompute_aux_precip!(::Precipitation0M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    #@. aux.precip_variables = (; q_rai = FT(0), q_sno = FT(0))
    aux.precip_variables.q_rai = FT(0)
    aux.precip_variables.q_sno = FT(0)

    tmp = @. precip_helper_sources_0M!(
        aux.params,
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
    microphys_params = KP.microphysics_params(aux.params)

    #@. aux.precip_variables = precip_helper_vars(Y.ρq_rai, Y.ρq_sno, aux.moisture_variables.ρ)
    #@. aux.precip_sources = precip_helper_sources_1M!(aux.params, aux.moisture_variables.ts, aux.moisture_variables.q_tot, aux.moisture_variables.q_liq, aux.moisture_variables.q_ice, aux.precip_variables.q_rai, aux.precip_variables.q_sno, aux.moisture_variables.T, aux.moisture_variables.ρ, aux.TS.dt)

    tmp_v = @. precip_helper_vars(Y.ρq_rai, Y.ρq_sno, aux.moisture_variables.ρ)
    aux.precip_variables.q_rai = tmp_v.q_rai
    aux.precip_variables.q_sno = tmp_v.q_sno

    tmp = @. precip_helper_sources_1M!(
        aux.params,
        aux.moisture_variables.ts,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.moisture_variables.q_ice,
        aux.precip_variables.q_rai,
        aux.precip_variables.q_sno,
        aux.moisture_variables.T,
        aux.moisture_variables.ρ,
        aux.TS.dt,
        ps,
    )
    aux.precip_sources.S_q_rai = tmp.S_q_rai
    aux.precip_sources.S_q_sno = tmp.S_q_sno
    aux.precip_sources.S_q_tot = tmp.S_q_tot
    aux.precip_sources.S_q_liq = tmp.S_q_liq
    aux.precip_sources.S_q_ice = tmp.S_q_ice

    If = CC.Operators.InterpolateC2F(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. aux.precip_velocities.term_vel_rai = CC.Geometry.WVector(If(tmp.term_vel_rai) * FT(-1))
    @. aux.precip_velocities.term_vel_sno = CC.Geometry.WVector(If(tmp.term_vel_sno) * FT(-1))
end
@inline function precompute_aux_precip!(ps::Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_rai)

    aux.aerosol_variables.N_aer = Y.N_aer
    tmp = @. aerosol_activation_helper(
        aux.params,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.aerosol_variables.N_aer,
        aux.aerosol_variables.N_aer_0,
        aux.moisture_variables.T,
        aux.moisture_variables.p,
        aux.moisture_variables.ρ,
        ClimaCore.Operators.InterpolateF2C().(aux.prescribed_velocity.ρw.components.data.:1),
        aux.TS.dt,
    )
    aux.aerosol_variables.S_N_aer = tmp.S_Na
    aux.precip_sources.S_N_liq = tmp.S_Nl

    aux.precip_variables.q_rai = Y.ρq_rai ./ aux.moisture_variables.ρ
    aux.precip_variables.N_liq = Y.N_liq
    aux.precip_variables.N_rai = Y.N_rai
    tmp = @. precip_helper_sources_2M!(
        aux.params,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.precip_variables.q_rai,
        aux.precip_variables.N_liq,
        aux.precip_variables.N_rai,
        aux.moisture_variables.T,
        aux.moisture_variables.ρ,
        aux.TS.dt,
        ps,
    )
    aux.precip_sources.S_q_rai = tmp.S_q_rai
    aux.precip_sources.S_q_tot = tmp.S_q_tot
    aux.precip_sources.S_q_liq = tmp.S_q_liq
    @. aux.precip_sources.S_N_liq += tmp.S_N_liq
    aux.precip_sources.S_N_rai = tmp.S_N_rai

    If = CC.Operators.InterpolateC2F(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. aux.precip_velocities.term_vel_N_rai = CC.Geometry.WVector(If(tmp.term_vel_N_rai) * FT(-1))
    @. aux.precip_velocities.term_vel_rai = CC.Geometry.WVector(If(tmp.term_vel_rai) * FT(-1))
end

"""
   Advection Equation: ∂ϕ/dt = -∂(vΦ)
"""
@inline function advection_tendency!(sm::AbstractMoistureStyle, dY, Y, aux, t)
    error("advection_tendency not implemented for a given $sm")
end
@inline function advection_tendency!(sp::AbstractPrecipitationStyle, dY, Y, aux, t)
    error("advection_tendency not implemented for a given $sp")
end
@inline function advection_tendency!(::EquilibriumMoisture, dY, Y, aux, t)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * aux.q_surf)),
        top = CC.Operators.Extrapolate(),
    )
    @. dY.ρq_tot += -∂(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) * If(Y.ρq_tot))

    if Bool(aux.params.qtot_flux_correction)
        fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
        @. dY.ρq_tot += fcc(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ), Y.ρq_tot)
    end

    return dY
end
@inline function advection_tendency!(::NonEquilibriumMoisture, dY, Y, aux, t)
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
    if Bool(aux.params.qtot_flux_correction)
        @. dY.ρq_tot += fcc(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ), Y.ρq_tot)
    end
    @. dY.ρq_liq += fcc(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ), Y.ρq_liq)
    @. dY.ρq_ice += fcc(aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ), Y.ρq_ice)


    return dY
end
@inline function advection_tendency!(::Union{NoPrecipitation, Precipitation0M}, dY, Y, aux, t) end
@inline function advection_tendency!(::Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())

    @. dY.ρq_rai +=
        -∂(
            (aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) + aux.precip_velocities.term_vel_rai) *
            If(Y.ρq_rai),
        )
    @. dY.ρq_sno +=
        -∂(
            (aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) + aux.precip_velocities.term_vel_sno) *
            If(Y.ρq_sno),
        )

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.ρq_rai +=
        fcc((aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) + aux.precip_velocities.term_vel_rai), Y.ρq_rai)
    @. dY.ρq_sno +=
        fcc((aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) + aux.precip_velocities.term_vel_sno), Y.ρq_sno)


    return dY
end
@inline function advection_tendency!(::Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())

    @. dY.N_rai +=
        -∂(
            (aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) + aux.precip_velocities.term_vel_N_rai) *
            If(Y.N_rai),
        )
    @. dY.ρq_rai +=
        -∂(
            (aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) + aux.precip_velocities.term_vel_rai) *
            If(Y.ρq_rai),
        )

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.N_rai +=
        fcc((aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) + aux.precip_velocities.term_vel_N_rai), Y.N_rai)
    @. dY.ρq_rai +=
        fcc((aux.prescribed_velocity.ρw / If(aux.moisture_variables.ρ) + aux.precip_velocities.term_vel_rai), Y.ρq_rai)

    return dY
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
    @. dY.N_liq += aux.precip_sources.S_N_liq
    @. dY.N_rai += aux.precip_sources.S_N_rai
    @. dY.N_aer += aux.aerosol_variables.S_N_aer

    return dY
end
