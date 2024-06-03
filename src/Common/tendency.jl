"""
   Different components of the ODE `rhs!` function,
   depending on the moisture and precipitation types.
"""
# some aliases to improve readability
const PP = TD.PhasePartition
const Lf = TD.latent_heat_fusion
const q_vap = TD.vapor_specific_humidity

# helper function to limit the tendency
function limit(q::FT, dt::FT, S::FT, n = 1) where {FT}
    return min(q / dt / n, S)
end

# helper function to safely get q from state ρq
function q_(ρq::FT, ρ::FT) where {FT}
    return max(FT(0), ρq / ρ)
end
@inline function zero_tendencies!(::CloudyPrecip, dY, Y, aux, t)
    FT = eltype(Y.ρq_vap)
    @. dY.N_aer = FT(0)
    @. dY.moments = ntuple(_ -> FT(0), length(dY.moments))
end

"""
    Zero out previous timestep tendencies
"""
@inline function zero_tendencies!(dY)
    @. dY = 0
end

"""
     Precompute the auxiliary values
"""
@inline function precompute_aux_thermo!(sm::AbstractMoistureStyle, Y, aux)
    error("precompute_aux not implemented for a given $sm")
end
@inline function precompute_aux_thermo!(::EquilibriumMoisture, Y, aux)

    (; thermo_params) = aux
    (; ts, ρ, ρ_dry, p, T, θ_dry, θ_liq_ice) = aux.thermo_variables
    (; q_tot, q_liq, q_ice) = aux.microph_variables

    @. ρ = ρ_dry + Y.ρq_tot
    @. ts = TD.PhaseEquil_ρTq(thermo_params, ρ, T, Y.ρq_tot / ρ)
    @. p = TD.air_pressure(thermo_params, ts)

    @. q_tot = TD.total_specific_humidity(thermo_params, ts)
    @. q_liq = TD.liquid_specific_humidity(thermo_params, ts)
    @. q_ice = TD.ice_specific_humidity(thermo_params, ts)

    @. θ_dry = TD.dry_pottemp(thermo_params, T, ρ_dry)
    @. θ_liq_ice = TD.liquid_ice_pottemp(thermo_params, ts)
end
@inline function precompute_aux_thermo!(::NonEquilibriumMoisture, Y, aux)

    (; thermo_params) = aux
    (; ts, ρ, ρ_dry, p, T, θ_dry, θ_liq_ice) = aux.thermo_variables
    (; q_tot, q_liq, q_ice) = aux.microph_variables

    @. ρ = ρ_dry + Y.ρq_tot

    @. q_tot = q_(Y.ρq_tot, ρ)
    @. q_liq = q_(Y.ρq_liq, ρ)
    @. q_ice = q_(Y.ρq_ice, ρ)

    @. ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, PP(q_tot, q_liq, q_ice))
    @. p = TD.air_pressure(thermo_params, ts)
    @. θ_dry = TD.dry_pottemp(thermo_params, T, ρ_dry)
    @. θ_liq_ice = TD.liquid_ice_pottemp(thermo_params, ts)
end
@inline function precompute_aux_thermo!(::CloudyMoisture, Y, aux)

    (; thermo_params) = aux
    (; ts, ρ, ρ_dry, p, T, θ_dry, θ_liq_ice) = aux.thermo_variables
    (; q_tot, q_liq, q_ice) = aux.microph_variables

    FT = eltype(Y.ρq_vap)
    @. q_tot = q_(Y.moments.:2 + Y.ρq_vap, ρ)
    @. q_liq = q_(Y.moments.:2, ρ)
    @. q_ice = FT(0)

    @. ρ = ρ_dry + Y.moments.:2 + Y.ρq_vap
    @. ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, PP(q_tot, q_liq, q_ice))
    @. p = TD.air_pressure(thermo_params, ts)
    @. θ_liq_ice = TD.liquid_ice_pottemp(thermo_params, ts)
    @. θ_dry = TD.dry_pottemp(thermo_params, T, ρ_dry)
end

@inline function precompute_aux_precip!(::Union{NoPrecipitation, Precipitation0M}, Y, aux) end
@inline function precompute_aux_precip!(ps::Precipitation1M, Y, aux)

    FT = eltype(Y.ρq_rai)

    (; ρ) = aux.thermo_variables
    (; q_rai, q_sno) = aux.microph_variables
    (; term_vel_rai, term_vel_sno) = aux.velocities

    @. q_rai = q_(Y.ρq_rai, ρ)
    @. q_sno = q_(Y.ρq_sno, ρ)

    @. term_vel_rai = CM1.terminal_velocity(ps.rain, ps.sedimentation.rain, ρ, q_rai)
    @. term_vel_sno = CM1.terminal_velocity(ps.snow, ps.sedimentation.snow, ρ, q_sno)
end
@inline function precompute_aux_precip!(ps::Precipitation2M, Y, aux)

    FT = eltype(Y.ρq_rai)
    sb2006 = ps.rain_formation
    vel_scheme = ps.sedimentation
    (; ρ) = aux.thermo_variables
    (; q_rai, N_rai, N_liq) = aux.microph_variables
    (; term_vel_rai, term_vel_N_rai) = aux.velocities

    @. q_rai = q_(Y.ρq_rai, ρ)
    @. N_rai = max(FT(0), Y.N_rai)
    @. N_liq = max(FT(0), Y.N_liq)

    @. term_vel_N_rai = getindex(CM2.rain_terminal_velocity(sb2006, vel_scheme, q_rai, ρ, N_rai), 1)
    @. term_vel_rai = getindex(CM2.rain_terminal_velocity(sb2006, vel_scheme, q_rai, ρ, N_rai), 2)
end
@inline function precompute_aux_precip!(ps::PrecipitationP3, Y, aux)

    FT = eltype(Y.ρq_rai)
    (; ρ) = aux.thermo_variables
    (; q_rai, q_ice, q_rim, B_rim, N_rai, N_ice, N_liq) = aux.microph_variables

    @. q_rai = q_(Y.ρq_rai, ρ)
    @. q_ice = q_(Y.ρq_ice, ρ)
    @. q_rim = q_(Y.ρq_rim, ρ)
    @. B_rim = q_(Y.ρq_rim, ρ)
    @. N_rai = max(FT(0), Y.N_rai)
    @. N_liq = max(FT(0), Y.N_liq)
    @. N_ice = max(FT(0), Y.N_ice)

    # TODO...
end

# helper function for precomputing aux precip for cloudy
function get_dists_moments(moments, NProgMoms::NTuple{ND, Int}) where {ND}
    FT = eltype(moments)
    return ntuple(ND) do i
        ntuple(3) do j
            ind = i == 1 ? 0 : sum(NProgMoms[1:(i - 1)])
            if j <= NProgMoms[i]
                moments[ind + j]
            else
                FT(0)
            end
        end
    end
end
@inline function get_updated_pdists(moments, old_pdists, cloudy_params)
    mom_normed = moments ./ cloudy_params.mom_norms
    mom_i = get_dists_moments(mom_normed, cloudy_params.NProgMoms)
    ntuple(length(old_pdists)) do i
        if old_pdists[i] isa CL.ParticleDistributions.GammaPrimitiveParticleDistribution
            CL.ParticleDistributions.update_dist_from_moments(
                old_pdists[i],
                mom_i[i],
                param_range = (; :k => (0.1, 10.0)),
            )
        else
            CL.ParticleDistributions.update_dist_from_moments(old_pdists[i], mom_i[i])
        end
    end
end
@inline function get_weighted_vt(moments, pdists, cloudy_params)
    FT = eltype(pdists[1])
    sed_flux = CL.Sedimentation.get_sedimentation_flux(pdists, cloudy_params.vel)
    weighted_vt = ntuple(length(moments)) do i
        ifelse(moments[i] > FT(0), -1 * sed_flux[i] * cloudy_params.mom_norms[i] / moments[i], FT(0))
    end
    return weighted_vt
end
@inline function precompute_aux_precip!(ps::CloudyPrecip, Y, aux)

    (; cloudy_params) = aux
    (; ρ) = aux.thermo_variables
    (; q_rai, N_rai, N_liq, pdists, moments) = aux.microph_variables
    (; weighted_vt) = aux.velocities
    (; nm_cloud) = aux.cloudy_variables

    _valof(::Val{NM}) where {NM} = NM
    NM1 = _valof(nm_cloud)
    rain_number_ind = NM1 + 1
    rain_mass_ind = NM1 + 2
    @. N_liq = Y.moments.:1
    @. N_rai = Y.moments.:($$rain_number_ind)
    @. q_rai = q_(Y.moments.:($$rain_mass_ind), ρ)
    @. moments = Y.moments

    @. pdists = get_updated_pdists(moments, pdists, cloudy_params)
    @. weighted_vt = get_weighted_vt(moments, pdists, cloudy_params)
end

@inline function precompute_aux_moisture_sources!(sm::AbstractMoistureStyle, aux)
    error("precompute_aux not implemented for a given $sm")
end
@inline function precompute_aux_moisture_sources!(::EquilibriumMoisture, aux) end
@inline function precompute_aux_moisture_sources!(ne::NonEquilibriumMoisture, aux)

    (; thermo_params) = aux
    (; ρ, T) = aux.thermo_variables
    (; q_tot, q_liq, q_ice) = aux.microph_variables

    S_q_liq = aux.scratch.tmp
    S_q_ice = aux.scratch.tmp2

    # helper type and wrapper to populate the tuple with sources
    S_eltype = eltype(aux.cloud_sources)
    to_sources(args...) = S_eltype(tuple(args...))

    @. S_q_liq = CMNe.conv_q_vap_to_q_liq_ice(
        ne.liquid,
        PP(thermo_params, TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)),
        PP(q_tot, q_liq, q_ice),
    )
    @. S_q_ice = CMNe.conv_q_vap_to_q_liq_ice(
        ne.ice,
        PP(thermo_params, TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)),
        PP(q_tot, q_liq, q_ice),
    )

    @. aux.cloud_sources = to_sources(S_q_liq, S_q_ice)
end
@inline function precompute_aux_moisture_sources!(::CloudyMoisture, dY, Y, aux, t) end

@inline function precompute_aux_precip_sources!(sp::AbstractPrecipitationStyle, aux)
    error("precompute_aux not implemented for a given $sp")
end
@inline function precompute_aux_precip_sources!(::NoPrecipitation, aux) end
@inline function precompute_aux_precip_sources!(ps::Precipitation0M, aux)

    (; thermo_params) = aux
    (; dt) = aux.TS
    (; ts) = aux.thermo_variables
    (; q_tot, q_liq, q_ice) = aux.microph_variables

    S_qt = aux.scratch.tmp
    λ = aux.scratch.tmp2

    @. λ = TD.liquid_fraction(thermo_params, ts)

    @. S_qt =
        -limit(
            q_liq + q_ice,
            dt,
            -CM0.remove_precipitation(ps.params, PP(q_tot, q_liq, q_ice), TD.q_vap_saturation(thermo_params, ts)),
        )

    S_eltype = eltype(aux.precip_sources)
    to_sources(args...) = S_eltype(tuple(args...))

    @. aux.precip_sources = to_sources(S_qt, S_qt * λ, S_qt * (1 - λ))
end
@inline function precompute_aux_precip_sources!(ps::Precipitation1M, aux)

    (; dt) = aux.TS
    (; thermo_params, common_params, air_params) = aux
    FT = eltype(thermo_params)
    (; ts, ρ, T) = aux.thermo_variables
    (; q_tot, q_liq, q_ice, q_rai, q_sno) = aux.microph_variables

    S = aux.scratch.tmp
    α = aux.scratch.tmp2

    T_fr::FT = TD.Parameters.T_freeze(thermo_params)
    c_vl::FT = TD.Parameters.cv_l(thermo_params)

    # helper type and wrapper to populate the tuple with sources
    S_eltype = eltype(aux.precip_sources)
    to_sources(args...) = S_eltype(tuple(args...))

    # zero out the aux sources
    @. aux.precip_sources = to_sources(0, 0, 0, 0, 0)
    #TODO - figure out a better way to split between equil and non-equil sources
    if Bool(common_params.precip_sources)
        rf = ps.rain_formation
        # autoconversion of liquid to rain
        if rf isa CMP.Acnv1M{FT}
            @. S = -limit(q_liq, dt, CM1.conv_q_liq_to_q_rai(rf, q_liq, true))
        elseif typeof(rf) in [CMP.LD2004{FT}, CMP.VarTimescaleAcnv{FT}]
            @. S = -limit(q_liq, dt, CM2.conv_q_liq_to_q_rai(rf, q_liq, ρ, common_params.prescribed_Nd))
        elseif typeof(rf.acnv) in [CMP.AcnvKK2000{FT}, CMP.AcnvB1994{FT}, CMP.AcnvTC1980{FT}]
            @. S = -limit(q_liq, dt, CM2.conv_q_liq_to_q_rai(rf, q_liq, ρ, common_params.prescribed_Nd))
        else
            error("Unrecognized rain formation scheme")
        end
        @. aux.precip_sources += to_sources(S, S, 0, -S, 0)

        # autoconversion of ice to snow
        @. S =
            -limit(q_ice, dt, CM1.conv_q_ice_to_q_sno(ps.ice, air_params, thermo_params, PP(q_tot, q_liq, q_ice), ρ, T))
        @. aux.precip_sources += to_sources(S, 0, S, 0, -S)

        # accretion cloud water + rain
        if typeof(rf) in [CMP.Acnv1M{FT}, CMP.LD2004{FT}, CMP.VarTimescaleAcnv{FT}]
            @. S = limit(q_liq, dt, CM1.accretion(ps.liquid, ps.rain, ps.sedimentation.rain, ps.ce, q_liq, q_rai, ρ))
        elseif typeof(rf.accr) in [CMP.AccrKK2000{FT}, CMP.AccrB1994{FT}]
            @. S = limit(q_liq, dt, CM2.accretion(rf, q_liq, q_rai, ρ))
        elseif rf.accr isa CMP.AccrTC1980{FT}
            @. S = limit(q_liq, dt, CM2.accretion(rf, q_liq, q_rai))
        else
            error("Unrecognized rain formation scheme")
        end
        @. aux.precip_sources += to_sources(-S, -S, 0, S, 0)

        # accretion cloud ice + snow
        @. S = limit(q_ice, dt, CM1.accretion(ps.ice, ps.snow, ps.sedimentation.snow, ps.ce, q_ice, q_sno, ρ))
        @. aux.precip_sources += to_sources(-S, 0, -S, 0, S)

        # sink of cloud water via accretion cloud water + snow
        @. S = -limit(q_liq, dt, CM1.accretion(ps.liquid, ps.snow, ps.sedimentation.snow, ps.ce, q_liq, q_sno, ρ))
        @. α = c_vl / Lf(thermo_params, ts) * (T - T_fr)
        @. aux.precip_sources += ifelse(T < T_fr, to_sources(S, S, 0, 0, -S), to_sources(S, S, 0, -S * (1 + α), S * α))

        # sink of cloud ice via accretion cloud ice - rain
        @. S = -limit(q_ice, dt, CM1.accretion(ps.ice, ps.rain, ps.sedimentation.rain, ps.ce, q_ice, q_rai, ρ))
        @. aux.precip_sources += to_sources(S, 0, S, 0, -S)

        # sink of rain via accretion cloud ice - rain
        @. S =
            -limit(q_rai, dt, CM1.accretion_rain_sink(ps.rain, ps.ice, ps.sedimentation.rain, ps.ce, q_ice, q_rai, ρ))
        @. aux.precip_sources += to_sources(0, 0, 0, S, -S)

        # accretion rain - snow
        @. S = ifelse(
            T < T_fr,
            limit(
                q_rai,
                dt,
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
            ),
            -limit(
                q_sno,
                dt,
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
            ),
        )
        @. aux.precip_sources += to_sources(0, 0, 0, -S, S)
    end

    if Bool(common_params.precip_sinks)
        # evaporation
        @. S =
            -limit(
                q_rai,
                dt,
                -CM1.evaporation_sublimation(
                    ps.rain,
                    ps.sedimentation.rain,
                    air_params,
                    thermo_params,
                    PP(q_tot, q_liq, q_ice),
                    q_rai,
                    ρ,
                    T,
                ),
            )
        @. aux.precip_sources += to_sources(-S, 0, 0, S, 0)

        # melting
        @. S = -limit(q_sno, dt, CM1.snow_melt(ps.snow, ps.sedimentation.snow, air_params, thermo_params, q_sno, ρ, T))
        @. aux.precip_sources += to_sources(0, 0, 0, -S, S)

        # deposition, sublimation
        @. S = CM1.evaporation_sublimation(
            ps.snow,
            ps.sedimentation.snow,
            air_params,
            thermo_params,
            PP(q_tot, q_liq, q_ice),
            q_sno,
            ρ,
            T,
        )
        @. S = ifelse(S > 0, limit(q_vap(thermo_params, ts), dt, S), -limit(q_sno, dt, -S))
        @. aux.precip_sources += to_sources(-S, 0, 0, 0, S)
    end
end
@inline function precompute_aux_precip_sources!(ps::Precipitation2M, aux)

    sb2006 = ps.rain_formation
    (; dt) = aux.TS
    (; thermo_params, common_params, air_params) = aux
    FT = eltype(thermo_params)
    (; ρ, T) = aux.thermo_variables
    (; q_tot, q_liq, q_ice, q_rai, N_liq, N_rai) = aux.microph_variables
    S₁ = aux.scratch.tmp
    S₂ = aux.scratch.tmp2

    # helper type and wrapper to populate the tuple with sources
    S_eltype = eltype(aux.precip_sources)
    to_sources(args...) = S_eltype(tuple(args...))

    # zero out the aux sources
    @. aux.precip_sources = to_sources(0, 0, 0, 0, 0, 0)
    #tot, liq, rai, Na, Nl, Nr
    if Bool(common_params.precip_sources)
        # autoconversion liquid to rain (mass)
        @. S₁ = limit(q_liq, dt, CM2.autoconversion(sb2006.acnv, q_liq, q_rai, ρ, N_liq).dq_rai_dt)
        @. aux.precip_sources += to_sources(-S₁, -S₁, S₁, 0, 0, 0)

        # autoconversion liquid to rain (number)
        @. S₂ = limit(N_liq, dt, CM2.autoconversion(sb2006.acnv, q_liq, q_rai, ρ, N_liq).dN_rai_dt, 2)
        @. aux.precip_sources += to_sources(0, 0, 0, 0, -2 * S₂, S₂)
        # liquid self_collection
        @. S₁ = -limit(N_liq, dt, -CM2.liquid_self_collection(sb2006.acnv, q_liq, ρ, -2 * S₂))
        @. aux.precip_sources += to_sources(0, 0, 0, 0, S₁, 0)

        # rain self_collection
        @. S₁ = CM2.rain_self_collection(sb2006.pdf_r, sb2006.self, q_rai, ρ, N_rai)
        @. aux.precip_sources += to_sources(0, 0, 0, 0, 0, -limit(N_rai, dt, -S₁))
        # rain breakup
        @. aux.precip_sources += to_sources(
            0,
            0,
            0,
            0,
            0,
            limit(N_rai, dt, CM2.rain_breakup(sb2006.pdf_r, sb2006.brek, q_rai, ρ, N_rai, S₁)),
        )

        # accretion cloud water + rain
        @. S₁ = limit(q_liq, dt, CM2.accretion(sb2006, q_liq, q_rai, ρ, N_liq).dq_rai_dt)
        @. S₂ = -limit(N_liq, dt, -CM2.accretion(sb2006, q_liq, q_rai, ρ, N_liq).dN_liq_dt)
        @. aux.precip_sources += to_sources(-S₁, -S₁, S₁, 0, S₂, 0)
    end

    if Bool(common_params.precip_sinks)
        # evaporation
        @. S₁ =
            -limit(
                N_rai,
                dt,
                -CM2.rain_evaporation(
                    sb2006,
                    air_params,
                    thermo_params,
                    PP(q_tot, q_liq, q_ice),
                    q_rai,
                    ρ,
                    N_rai,
                    T,
                ).evap_rate_0,
            )
        @. S₂ =
            -limit(
                q_rai,
                dt,
                -CM2.rain_evaporation(
                    sb2006,
                    air_params,
                    thermo_params,
                    PP(q_tot, q_liq, q_ice),
                    q_rai,
                    ρ,
                    N_rai,
                    T,
                ).evap_rate_1,
            )
        @. aux.precip_sources += to_sources(-S₂, 0, S₂, 0, 0, S₁)
    end
end

@inline function precompute_aux_precip_sources!(ps::PrecipitationP3, aux)
    return nothing
end

# helper function for precomputing aux sources for cloudy
@inline function get_coal_sources(cloudy_params, moments, pdists, dt)

    dY_coal_tmp = CL.Coalescence.get_coal_ints(CL.EquationTypes.AnalyticalCoalStyle(), pdists, cloudy_params.coal_data)
    S_coal = ntuple(length(moments)) do j
        ifelse(
            dY_coal_tmp[j] >= 0,
            dY_coal_tmp[j] * cloudy_params.mom_norms[j],
            max(dY_coal_tmp[j] * cloudy_params.mom_norms[j], -moments[j] / dt),
        )
    end

    return S_coal
end
@inline function get_cond_evap_sources(
    thermo_params,
    air_params,
    cloudy_params,
    q_tot,
    q_liq,
    moments,
    pdists,
    T,
    ρ,
    dt,
)
    FT = eltype(pdists[1])
    q = TD.PhasePartition(q_tot, q_liq, FT(0))

    ξ = CM.Common.G_func(air_params, thermo_params, T, TD.Liquid())
    ξ_normed = ξ / cloudy_params.norms[2]^(2 / 3)
    s = TD.supersaturation(thermo_params, q, ρ, T, TD.Liquid())
    dY_ce_tmp = CL.Condensation.get_cond_evap(pdists, s, ξ_normed) .* cloudy_params.mom_norms
    S_cond_evap = ntuple(length(moments)) do j
        if (j == cloudy_params.NProgMoms[1]+1 && moments[j+1] / moments[j] < 1e-11) || (j == 1 && moments[j+1] / moments[j] < 1e-14)
            -moments[j] / dt / 4
        else
            ifelse(dY_ce_tmp[j] >= 0, dY_ce_tmp[j], max(dY_ce_tmp[j], -moments[j] / dt))
        end
    end

    return S_cond_evap
end
@inline function precompute_aux_precip_sources!(ps::CloudyPrecip, aux)

    (; common_params, thermo_params, air_params, cloudy_params) = aux
    (; q_tot, q_liq, moments, pdists) = aux.microph_variables
    (; T, ρ) = aux.thermo_variables
    (; dt) = aux.TS
    (; tmp_cloudy) = aux.scratch
    (; nm_cloud) = aux.cloudy_variables

    FT = eltype(thermo_params)
    @. aux.precip_sources.moments *= FT(0)

    # coalescence
    if Bool(common_params.precip_sources)
        @. aux.precip_sources.moments += get_coal_sources(cloudy_params, moments, pdists, dt)
    end

    # condensation or evaporations
    if Bool(common_params.precip_sinks)
        @. tmp_cloudy =
            get_cond_evap_sources(thermo_params, air_params, cloudy_params, q_tot, q_liq, moments, pdists, T, ρ, dt)
        @. aux.precip_sources.moments += tmp_cloudy

        _valof(::Val{NM}) where {NM} = NM
        NM1 = _valof(nm_cloud)
        rain_mass_ind = NM1 + 2
        @. aux.precip_sources.ρq_vap = -tmp_cloudy.:2 - tmp_cloudy.:($$rain_mass_ind)
    end

end

"""
   Additional source terms
"""
@inline function cloud_sources_tendency!(ms::AbstractMoistureStyle, dY, Y, aux, t)
    error("sources_tendency not implemented for a given $ms")
end
@inline function cloud_sources_tendency!(::EquilibriumMoisture, dY, Y, aux, t) end
@inline function cloud_sources_tendency!(ms::NonEquilibriumMoisture, dY, Y, aux, t)

    precompute_aux_moisture_sources!(ms, aux)

    @. dY.ρq_liq += aux.thermo_variables.ρ * aux.cloud_sources.q_liq
    @. dY.ρq_ice += aux.thermo_variables.ρ * aux.cloud_sources.q_ice
    return dY
end
@inline function cloud_sources_tendency!(::MoistureP3, dY, Y, aux, t) end
@inline function cloud_sources_tendency!(::CloudyMoisture, dY, Y, aux, t)
    @. dY.ρq_vap += aux.precip_sources.ρq_vap
    @. dY.ρq_vap += aux.activation_sources.ρq_vap
    return dY
end


@inline function precip_sources_tendency!(ms::AbstractMoistureStyle, ps::AbstractPrecipitationStyle, dY, Y, aux, t)
    error("sources_tendency not implemented for a given $sp")
end
@inline function precip_sources_tendency!(ms::AbstractMoistureStyle, ::NoPrecipitation, dY, Y, aux, t) end

@inline function precip_sources_tendency!(ms::AbstractMoistureStyle, ps::Precipitation0M, dY, Y, aux, t)

    precompute_aux_precip_sources!(ps, aux)

    @. dY.ρq_tot += aux.thermo_variables.ρ * aux.precip_sources.q_tot

    if ms isa NonEquilibriumMoisture
        @. dY.ρq_liq += aux.thermo_variables.ρ * aux.precip_sources.q_liq
        @. dY.ρq_ice += aux.thermo_variables.ρ * aux.precip_sources.q_ice
    end
end

@inline function precip_sources_tendency!(ms::AbstractMoistureStyle, ps::Precipitation1M, dY, Y, aux, t)

    precompute_aux_precip_sources!(ps, aux)

    @. dY.ρq_tot += aux.thermo_variables.ρ * aux.precip_sources.q_tot
    @. dY.ρq_rai += aux.thermo_variables.ρ * aux.precip_sources.q_rai
    @. dY.ρq_sno += aux.thermo_variables.ρ * aux.precip_sources.q_sno

    if ms isa NonEquilibriumMoisture
        @. dY.ρq_liq += aux.thermo_variables.ρ * aux.precip_sources.q_liq
        @. dY.ρq_ice += aux.thermo_variables.ρ * aux.precip_sources.q_ice
    end

    return dY
end
@inline function precip_sources_tendency!(ms::AbstractMoistureStyle, ps::Precipitation2M, dY, Y, aux, t)

    precompute_aux_precip_sources!(ps, aux)

    @. dY.ρq_tot += aux.thermo_variables.ρ * aux.precip_sources.q_tot
    @. dY.ρq_rai += aux.thermo_variables.ρ * aux.precip_sources.q_rai
    @. dY.N_liq += aux.precip_sources.N_liq
    @. dY.N_rai += aux.precip_sources.N_rai

    if ms isa NonEquilibriumMoisture
        @. dY.ρq_liq += aux.thermo_variables.ρ * aux.precip_sources.q_liq
    end

    @. dY.N_liq += aux.activation_sources.N_liq
    @. dY.N_aer += aux.activation_sources.N_aer

    return dY
end
@inline function precip_sources_tendency!(ms::MoistureP3, ps::PrecipitationP3, dY, Y, aux, t)
    return dY
end
@inline function precip_sources_tendency!(ms::CloudyMoisture, ps::CloudyPrecip, dY, Y, aux, t)

    precompute_aux_precip_sources!(ps, aux)

    @. dY.moments += aux.precip_sources.moments
    @. dY.moments += aux.activation_sources.activation
    @. dY.N_aer += aux.activation_sources.N_aer
    return dY
end
