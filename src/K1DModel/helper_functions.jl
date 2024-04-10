"""
    Prescribed momentum flux as a function of time
"""
@inline function ρw_helper(t, w1, t1)
    return t < t1 ? w1 * sin(pi * t / t1) : 0.0
end

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

"""
    Returns mositure and precipitation types
"""
function get_moisture_and_precipitation_types(
    FT,
    moisture_choice::String,
    precipitation_choice::String,
    rain_formation_choice::String,
    sedimentation_choice::String,
    toml_dict,
)
    if moisture_choice == "EquilibriumMoisture"
        moisture = CO.EquilibriumMoisture_ρdTq()
    elseif moisture_choice == "NonEquilibriumMoisture"
        moisture = CO.NonEquilibriumMoisture_ρdTq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
    else
        error("Invalid moisture choice: $moisture_choice")
    end
    if precipitation_choice == "NoPrecipitation"
        precip = CO.NoPrecipitation()
    elseif precipitation_choice == "Precipitation0M"
        precip = CO.Precipitation0M(CMP.Parameters0M(FT, toml_dict))
    elseif precipitation_choice == "Precipitation1M"
        if sedimentation_choice == "CliMA_1M"
            st = CMP.Blk1MVelType(FT, toml_dict)
        elseif sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelType(FT, toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "CliMA_1M"
            rain_params = CMP.Rain(FT, toml_dict)
            rf = rain_params.acnv1M
        elseif rain_formation_choice == "KK2000"
            rf = CMP.KK2000(FT, toml_dict)
        elseif rain_formation_choice == "B1994"
            rf = CMP.B1994(FT, toml_dict)
        elseif rain_formation_choice == "TC1980"
            rf = CMP.TC1980(FT, toml_dict)
        elseif rain_formation_choice == "LD2004"
            rf = CMP.LD2004(FT, toml_dict)
        elseif rain_formation_choice == "VarTimeScaleAcnv"
            rf = CMP.VarTimescaleAcnv(FT, toml_dict)
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
        precip = CO.Precipitation1M(
            CMP.CloudLiquid(FT, toml_dict),
            CMP.CloudIce(FT, toml_dict),
            CMP.Rain(FT, toml_dict),
            CMP.Snow(FT, toml_dict),
            CMP.CollisionEff(FT, toml_dict),
            rf,
            st,
        )
    elseif precipitation_choice == "Precipitation2M"
        if sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelTypeRain(FT, toml_dict)
        elseif sedimentation_choice == "SB2006"
            st = CMP.SB2006VelType(FT, toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "SB2006"
            precip = CO.Precipitation2M(CMP.SB2006(FT, toml_dict), st)
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
    else
        error("Invalid precipitation choice: $precipitation_choice")
    end
    return moisture, precip
end

"""
    Returns the field of a given variable from ODE solutions
"""
function get_variable_data_from_ODE(u, aux, precip, var::String)

    ρ = parent(aux.moisture_variables.ρ_dry)[:] .+ parent(u.ρq_tot)[:]

    if var == "qt"
        output = parent(u.ρq_tot) ./ ρ
    elseif var == "ql"
        output = parent(u.ρq_liq) ./ ρ
    elseif var == "qr"
        output = parent(u.ρq_rai) ./ ρ
    elseif var == "qv"
        output = (parent(u.ρq_tot) .- parent(u.ρq_liq)) ./ ρ
    elseif var == "qlr"
        output = (parent(u.ρq_liq) .+ parent(u.ρq_rai)) ./ ρ
    elseif var == "rt"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = _qtot ./ (1 .- _qtot)
    elseif var == "rl"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = parent(u.ρq_liq) ./ ρ ./ (1 .- _qtot)
    elseif var == "rr"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = parent(u.ρq_rai) ./ ρ ./ (1 .- _qtot)
    elseif var == "rv"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = (parent(u.ρq_tot) .- parent(u.ρq_liq)) ./ ρ ./ (1 .- _qtot)
    elseif var == "rlr"
        _qtot = parent(u.ρq_tot) ./ ρ
        output = (parent(u.ρq_liq) .+ parent(u.ρq_rai)) ./ ρ ./ (1 .- _qtot)
    elseif var == "Nl"
        output = parent(u.N_liq)
    elseif var == "Nr"
        output = parent(u.N_rai)
    elseif var == "Na"
        output = parent(u.N_aer)
    elseif var == "rho"
        output = ρ
    elseif var == "rain averaged terminal velocity"
        qr = parent(u.ρq_rai) ./ ρ
        if precip isa Precipitation1M
            vt = CM1.terminal_velocity.(precip.rain, precip.sedimentation.rain, ρ, qr)
        elseif precip isa Precipitation2M
            Nr = parent(u.N_rai)
            f1(xx, yy, zz) = CM2.rain_terminal_velocity(precip.rain_formation, precip.sedimentation, xx, yy, zz)[2]
            vt = f1.(qr, ρ, Nr)
        else
            error("Computing rain averaged terminal velocity for the given precipitation style is invalid!!")
        end
        output = vt
    elseif var == "rainrate"
        qr = parent(u.ρq_rai) ./ ρ
        if precip isa Precipitation1M
            vt = CM1.terminal_velocity.(precip.rain, precip.sedimentation.rain, ρ, qr)
        elseif precip isa Precipitation2M
            Nr = parent(u.N_rai)
            f2(xx, yy, zz) = CM2.rain_terminal_velocity(precip.rain_formation, precip.sedimentation, xx, yy, zz)[2]
            vt = f2.(qr, ρ, Nr)
        else
            error("Computing rainrate for the given precipitation style is invalid!!")
        end
        output = qr .* ρ .* vt .* 3600
    else
        error("Data name \"" * var * "\" not recognized!!")
    end

    return output

end
