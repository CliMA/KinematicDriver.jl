"""
    Returns moisture type
"""
function get_moisture_type(
    FT,
    moisture_choice::String,
    toml_dict,
)
    if moisture_choice == "EquilibriumMoisture"
        moisture = EquilibriumMoisture_ρdTq()
    elseif moisture_choice == "NonEquilibriumMoisture"
        moisture = NonEquilibriumMoisture_ρdTq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
    else
        error("Invalid moisture choice: $moisture_choice")
    end
    return moisture
end

"""
    Returns precipitation type
"""
function get_precipitation_type(
    FT,
    precipitation_choice::String,
    toml_dict;
    rain_formation_choice::Union{Nothing, String} = nothing,
    sedimentation_choice::Union{Nothing, String} = nothing
)
    if precipitation_choice == "NoPrecipitation"
        precip = NoPrecipitation()
    elseif precipitation_choice == "Precipitation0M"
        precip = Precipitation0M(CMP.Parameters0M(FT, toml_dict))
    elseif precipitation_choice == "Precipitation1M"
        if sedimentation_choice == "CliMA_1M" || sedimentation_choice === nothing
            st = CMP.Blk1MVelType(FT, toml_dict)
        elseif sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelType(FT, toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "CliMA_1M" || rain_formation_choice === nothing
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
        precip = Precipitation1M(
            CMP.CloudLiquid(FT, toml_dict),
            CMP.CloudIce(FT, toml_dict),
            CMP.Rain(FT, toml_dict),
            CMP.Snow(FT, toml_dict),
            CMP.CollisionEff(FT, toml_dict),
            rf,
            st,
        )
    elseif precipitation_choice == "Precipitation2M"
        if sedimentation_choice == "SB2006" || sedimentation_choice === nothing
            st = CMP.SB2006VelType(FT, toml_dict)
        elseif sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelTypeRain(FT, toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "SB2006" || rain_formation_choice === nothing
            precip = Precipitation2M(CMP.SB2006(FT, toml_dict), st)
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
    else
        error("Invalid precipitation choice: $precipitation_choice")
    end
    return precip
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