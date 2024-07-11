"""
    Returns moisture type
"""
function get_moisture_type(moisture_choice::String, toml_dict)
    if moisture_choice == "EquilibriumMoisture"
        moisture = EquilibriumMoisture()
    elseif moisture_choice == "NonEquilibriumMoisture"
        moisture = NonEquilibriumMoisture(CMP.CloudLiquid(toml_dict), CMP.CloudIce(toml_dict))
    elseif moisture_choice == "CloudyMoisture"
        moisture = CloudyMoisture()
    else
        error("Invalid moisture choice: $moisture_choice")
    end
    return moisture
end

"""
    Returns precipitation type
"""
function get_precipitation_type(
    precipitation_choice::String,
    toml_dict;
    rain_formation_choice::Union{Nothing, String} = nothing,
    sedimentation_choice::Union{Nothing, String} = nothing,
)
    if precipitation_choice == "NoPrecipitation"
        precip = NoPrecipitation()
    elseif precipitation_choice == "Precipitation0M"
        precip = Precipitation0M(CMP.Parameters0M(toml_dict))
    elseif precipitation_choice == "Precipitation1M"
        if sedimentation_choice == "CliMA_1M" || sedimentation_choice === nothing
            st = CMP.Blk1MVelType(toml_dict)
        elseif sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelType(toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "CliMA_1M" || rain_formation_choice === nothing
            rain_params = CMP.Rain(toml_dict)
            rf = rain_params.acnv1M
        elseif rain_formation_choice == "KK2000"
            rf = CMP.KK2000(toml_dict)
        elseif rain_formation_choice == "B1994"
            rf = CMP.B1994(toml_dict)
        elseif rain_formation_choice == "TC1980"
            rf = CMP.TC1980(toml_dict)
        elseif rain_formation_choice == "LD2004"
            rf = CMP.LD2004(toml_dict)
        elseif rain_formation_choice == "VarTimeScaleAcnv"
            rf = CMP.VarTimescaleAcnv(toml_dict)
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
        precip = Precipitation1M(
            CMP.CloudLiquid(toml_dict),
            CMP.CloudIce(toml_dict),
            CMP.Rain(toml_dict),
            CMP.Snow(toml_dict),
            CMP.CollisionEff(toml_dict),
            rf,
            st,
        )
    elseif precipitation_choice == "Precipitation2M"
        if sedimentation_choice == "SB2006" || sedimentation_choice === nothing
            st = CMP.SB2006VelType(toml_dict)
        elseif sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelTypeRain(toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "SB2006" || rain_formation_choice === nothing
            precip = Precipitation2M(CMP.SB2006(toml_dict, false), st)
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
    elseif precipitation_choice == "CloudyPrecip"
        precip = CloudyPrecip()
    else
        error("Invalid precipitation choice: $precipitation_choice")
    end
    return precip
end

"""
    Returns the field of a given variable from ODE solutions
"""
function get_variable_data_from_ODE(u, aux, precip, var::String)

    ρ = parent(aux.thermo_variables.ρ_dry)[:] .+ parent(u.ρq_tot)[:]

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
        _qrai = parent(u.ρq_rai) ./ ρ
        if precip isa Precipitation1M
            vt = CM1.terminal_velocity.(precip.rain, precip.sedimentation.rain, ρ, _qrai)
        elseif precip isa Precipitation2M
            _Nrai = parent(u.N_rai)
            f1(xx, yy, zz) = CM2.rain_terminal_velocity(precip.rain_formation, precip.sedimentation, xx, yy, zz)[2]
            vt = f1.(_qrai, ρ, _Nrai)
        else
            error("Computing rain averaged terminal velocity for the given precipitation style is invalid!!")
        end
        output = vt
    elseif var == "rainrate"
        _qrai = parent(u.ρq_rai) ./ ρ
        if precip isa Precipitation1M
            vt = CM1.terminal_velocity.(precip.rain, precip.sedimentation.rain, ρ, _qrai)
        elseif precip isa Precipitation2M
            _Nrai = parent(u.N_rai)
            f2(xx, yy, zz) = CM2.rain_terminal_velocity(precip.rain_formation, precip.sedimentation, xx, yy, zz)[2]
            vt = f2.(_qrai, ρ, _Nrai)
        else
            error("Computing rainrate for the given precipitation style is invalid!!")
        end
        output = _qrai .* ρ .* vt .* 3600
    elseif var == "Z"
        _qliq = parent(u.ρq_liq) ./ ρ
        _qrai = parent(u.ρq_rai) ./ ρ
        _Nliq = parent(u.N_liq)
        _Nrai = parent(u.N_rai)

        if precip isa Precipitation1M
            output = CM1.radar_reflectivity.(precip.rain, _qrai, ρ)
        elseif precip isa Precipitation2M
            output = CM2.radar_reflectivity.(precip.rain_formation, _qliq, _qrai, _Nliq, _Nrai, ρ)
        else
            error("Computing radar reflectivity for the given precipitation style is invalid!!")
        end
    elseif var == "reff"
        _qliq = parent(u.ρq_liq) ./ ρ
        _qrai = parent(u.ρq_rai) ./ ρ
        _Nliq = parent(u.N_liq)
        _Nrai = parent(u.N_rai)

        if precip isa Precipitation2M
            output = CM2.effective_radius.(precip.rain_formation, _qliq, _qrai, _Nliq, _Nrai, ρ)
        else
            error("Computing effective radius for the given precipitation style is invalid!!")
        end
    elseif var == "rainrate_surface"
        _qrai = parent(u.ρq_rai) ./ ρ
        if precip isa Precipitation1M
            vt = CM1.terminal_velocity.(precip.rain, precip.sedimentation.rain, ρ, _qrai)
        elseif precip isa Precipitation2M
            _Nrai = parent(u.N_rai)
            f(xx, yy, zz) = CM2.rain_terminal_velocity(precip.rain_formation, precip.sedimentation, xx, yy, zz)[2]
            vt = f.(_qrai, ρ, _Nrai)
        else
            error("Computing rainrate for the given precipitation style is invalid!!")
        end
        rainrate_ = _qrai .* ρ .* vt .* 3600
        rainrate = 1.5 .* rainrate_[1] - 0.5 .* rainrate_[2]
        output = [max(0.0, rainrate)]
    elseif var == "reff_top"
        _qliq = parent(u.ρq_liq) ./ ρ
        _qrai = parent(u.ρq_rai) ./ ρ
        _Nliq = parent(u.N_liq)
        _Nrai = parent(u.N_rai)

        z_top = length(_qliq)
        threshold = 1e-6
        while _qliq[z_top] < threshold && z_top > 1
            z_top -= 1 
        end     

        depth = 500.0
        z = parent(CC.Fields.coordinate_field(aux.prescribed_velocity.ρw).z)
        dz = (maximum(z) - minimum(z)) / length(z)
        steps = Int(round(depth / dz))
        _ind = max(z_top - steps + 1, 1)

        ql = mean(_qliq[_ind:z_top])
        qr = mean(_qrai[_ind:z_top])
        Nl = mean(_Nliq[_ind:z_top])
        Nr = mean(_Nrai[_ind:z_top])
        ρ = mean(ρ[_ind:z_top])
        if precip isa Precipitation2M
                output = [CM2.effective_radius(precip.rain_formation, ql, qr, Nl, Nr, ρ)]
        else
            error("Computing effectve radius for the given precipitation style is invalid!!")
        end
    elseif var == "Z_top"
        _qliq = parent(u.ρq_liq) ./ ρ
        _qrai = parent(u.ρq_rai) ./ ρ
        _Nliq = parent(u.N_liq)
        _Nrai = parent(u.N_rai)

        z_top = length(_qliq)
        threshold = 1e-6
        while _qliq[z_top] < threshold && z_top > 1
            z_top -= 1
        end

        depth = 500.0
        z = parent(CC.Fields.coordinate_field(aux.prescribed_velocity.ρw).z)
        dz = (maximum(z) - minimum(z)) / length(z)
        steps = Int(round(depth / dz))
        _ind = max(z_top - steps + 1, 1)

        ql = mean(_qliq[_ind:z_top])
        qr = mean(_qrai[_ind:z_top])
        Nl = mean(_Nliq[_ind:z_top])
        Nr = mean(_Nrai[_ind:z_top])
        ρ = mean(ρ[_ind:z_top])
        if precip isa Precipitation1M
            output = [CM1.radar_reflectivity(precip.rain, qr, ρ)]
        elseif precip isa Precipitation2M
            output = [CM2.radar_reflectivity(precip.rain_formation, ql, qr, Nl, Nr, ρ)]
            #output = output < -50.0 ? -50.0 : output
        else
            error("Computing radar reflectivity for the given precipitation style is invalid!!")
        end
    elseif var == "Z_mid"
        _qliq = parent(u.ρq_liq) ./ ρ
        _qrai = parent(u.ρq_rai) ./ ρ
        _Nliq = parent(u.N_liq)
        _Nrai = parent(u.N_rai)

        z_top = length(_qliq)
        threshold = 1e-6
        while _qliq[z_top] < threshold && z_top > 1
            z_top -= 1
        end
        
        depth = 500.0
        height = 500.0
        z = parent(CC.Fields.coordinate_field(aux.prescribed_velocity.ρw).z)
        dz = (maximum(z) - minimum(z)) / length(z)
        steps = Int(round(depth / (2 * dz)))
        _bottom_ind = Int(round(height / dz))
        _ind = Int(round((z_top + _bottom_ind) / 2))
        _ind1 = max(_ind - steps + 1, 1)
        _ind2 = min(_ind + steps - 1, z_top)

        ql = mean(_qliq[_ind1:_ind2])
        qr = mean(_qrai[_ind1:_ind2])
        Nl = mean(_Nliq[_ind1:_ind2])
        Nr = mean(_Nrai[_ind1:_ind2])
        ρ = mean(ρ[_ind1:_ind2])
        if precip isa Precipitation1M
            output = [CM1.radar_reflectivity(precip.rain, qr, ρ)]
        elseif precip isa Precipitation2M
            output = [CM2.radar_reflectivity(precip.rain_formation, ql, qr, Nl, Nr, ρ)]
            #output = output < -50.0 ? -50.0 : output
        else
            error("Computing radar reflectivity for the given precipitation style is invalid!!")
        end
    elseif var == "Z_bottom"
        _qliq = parent(u.ρq_liq) ./ ρ
        _qrai = parent(u.ρq_rai) ./ ρ
        _Nliq = parent(u.N_liq)
        _Nrai = parent(u.N_rai)

        z_top = length(_qliq)
        threshold = 1e-6
        while _qliq[z_top] < threshold && z_top > 1
            z_top -= 1
        end

        depth = 500.0
        height = 500.0
        z = parent(CC.Fields.coordinate_field(aux.prescribed_velocity.ρw).z)
        dz = (maximum(z) - minimum(z)) / length(z)
        steps = Int(round(depth / dz))
        _ind = Int(round(height / dz))
        _final_ind = min(_ind + steps - 1, z_top)

        ql = mean(_qliq[_ind:_final_ind])
        qr = mean(_qrai[_ind:_final_ind])
        Nl = mean(_Nliq[_ind:_final_ind])
        Nr = mean(_Nrai[_ind:_final_ind])
        ρ = mean(ρ[_ind:_final_ind])
        if precip isa Precipitation1M
            output = [CM1.radar_reflectivity(precip.rain, qr, ρ)]
        elseif precip isa Precipitation2M
            output = [CM2.radar_reflectivity(precip.rain_formation, ql, qr, Nl, Nr, ρ)]
            #output = output < -50.0 ? -50.0 : output
        else
            error("Computing radar reflectivity for the given precipitation style is invalid!!")
        end
    elseif var == "Z_1"
        _qliq = parent(u.ρq_liq) ./ ρ
        _qrai = parent(u.ρq_rai) ./ ρ
        _Nliq = parent(u.N_liq)
        _Nrai = parent(u.N_rai)

        z_top = length(_qliq)
        threshold = 1e-6
        while _qliq[z_top] < threshold && z_top > 1
            z_top -= 1
        end 

        depth = 500.0
        z = parent(CC.Fields.coordinate_field(aux.prescribed_velocity.ρw).z)
        dz = (maximum(z) - minimum(z)) / length(z)
        steps = Int(round(depth / dz))
        z_1 =  z_top - steps
        _ind = max(z_1 - steps + 1, 1)

        ql = mean(_qliq[_ind:z_1])
        qr = mean(_qrai[_ind:z_1])
        Nl = mean(_Nliq[_ind:z_1])
        Nr = mean(_Nrai[_ind:z_1])
        ρ = mean(ρ[_ind:z_1])
        if precip isa Precipitation1M
            output = [CM1.radar_reflectivity(precip.rain, qr, ρ)]
        elseif precip isa Precipitation2M
            output = [CM2.radar_reflectivity(precip.rain_formation, ql, qr, Nl, Nr, ρ)]
            #output = output < FT(-50) ? FT(-50) : output
        else
            error("Computing radar reflectivity for the given precipitation style is invalid!!")
        end
    elseif var == "Z_2"
        _qliq = parent(u.ρq_liq) ./ ρ
        _qrai = parent(u.ρq_rai) ./ ρ
        _Nliq = parent(u.N_liq)
        _Nrai = parent(u.N_rai)

        z_top = length(_qliq)
        threshold = 1e-6
        while _qliq[z_top] < threshold && z_top > 1
            z_top -= 1
        end

        depth = 500.0
        z = parent(CC.Fields.coordinate_field(aux.prescribed_velocity.ρw).z)
        dz = (maximum(z) - minimum(z)) / length(z)
        steps = Int(round(depth / dz))
        z_2 = z_top - (2 * steps)
        _ind = max(z_2 - steps + 1, 1)

        ql = mean(_qliq[_ind:z_2])
        #println(ql)
        qr = mean(_qrai[_ind:z_2])
        #println(qr)
        Nl = mean(_Nliq[_ind:z_2])
        #println(Nl)
        Nr = mean(_Nrai[_ind:z_2])
        #println(Nr)
        ρ = mean(ρ[_ind:z_2])
        #println(ρ)
        if precip isa Precipitation1M
            output = [CM1.radar_reflectivity(precip.rain, qr, ρ)]
            #println(output)
        elseif precip isa Precipitation2M
            output = [CM2.radar_reflectivity(precip.rain_formation, ql, qr, Nl, Nr, ρ)]
            #output = output < FT(-50) ? FT(-50) : output
            #println("Z_2")
            #println(output)
        else
            error("Computing radar reflectivity for the given precipitation style is invalid!!")
        end
    else
        error("Data name \"" * var * "\" not recognized!!")
    end

    return output

end
