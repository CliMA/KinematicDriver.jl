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
    elseif moisture_choice == "MoistureP3"
        moisture = MoistureP3()
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
    boundary::Union{Nothing, NamedTuple} = nothing,
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
            precip = Precipitation2M(CMP.SB2006(toml_dict), st)
        elseif rain_formation_choice == "SB2006NL"
            precip = Precipitation2M(CMP.SB2006(toml_dict, false), st)
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
    elseif precipitation_choice == "CloudyPrecip"
        precip = CloudyPrecip()
    elseif precipitation_choice == "PrecipitationP3"
        FT = Float64
        p3_params = CMP.ParametersP3(FT)
        Chen2022 = CMP.Chen2022VelType(FT)
        sb2006 = CMP.SB2006(toml_dict)
        precip = PrecipitationP3(p3_params, Chen2022, sb2006, boundary)
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
    elseif var == "Z_top"
        _qliq = parent(u.ρq_liq) ./ ρ
        _qrai = parent(u.ρq_rai) ./ ρ
        _qlr = _qliq + _qrai

        z_top = length(_qliq)
        threshold = 1e-3
        while sum(_qlr[z_top, :]) < threshold && z_top > 1
            z_top -= 1
        end

        depth = 500.0
        z = parent(CC.Fields.coordinate_field(aux.prescribed_velocity.ρw).z)
        dz = (maximum(z) - minimum(z)) / length(z)
        steps = Int(round(depth / dz))
        _ind = max(z_top - steps + 1, 1)

        qr = vec(Statistics.mean(_qrai[_ind:z_top, :], dims = 1))
        ρ = vec(Statistics.mean(ρ[_ind:z_top, :], dims = 1))

        if precip isa Precipitation1M
            output = CMD.radar_reflectivity_1M.(precip.rain, qr, ρ)
            output = replace!(output, -Inf => -300.0, Inf => 300.0, NaN => 0.0) # TODO this should be captured in CM1
        else
            error("Computing radar reflectivity for the given precipitation style is invalid!!")
        end
    else
        error("Data name \"" * var * "\" not recognized!!")
    end

    return output

end

"""
    triangle_inequality_limiter(force, limit)

Limits a `force` (or source term) based on a maximum allowable `limit` using a
formula derived from the triangle inequality, as proposed by Horn (2012).

The formula used is: `L = F + M - sqrt(F² + M²)`, where `F` is the `force` and
`M` is the `limit`.

This limiter is designed to smoothly reduce the `force` as it approaches or
exceeds the `limit`, ensuring the result `L` satisfies `0 ≤ L ≤ M` if `F ≥ 0`
and `M > 0`. It also preserves `L ≤ F`. It's particularly useful for ensuring
that source terms (e.g., emissions, chemical production rates) do not become
unphysically large or lead to numerical instability, while being continuously
differentiable.

Arguments:
- `force`: The original force or source term value.
- `limit`: The maximum permissible positive value for the limited force.

Returns:
- The limited force value.

Reference:
- Horn, M. (2012). "ASAMgpu V1.0 – a moist fully compressible atmospheric model using 
    graphics processing units (GPUs)". Geoscientific Model Development,
    5, 345–353. https://doi.org/10.5194/gmd-5-345-2012
"""

function triangle_inequality_limiter(force, limit)
    return force + limit - sqrt(force^2 + limit^2)
end
