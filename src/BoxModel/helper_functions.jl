"""
    Returns mositure and precipitation types
"""
function get_precipitation_type(FT, precipitation_choice::String, rain_formation_choice::String, toml_dict)
    if precipitation_choice == "Precipitation1M"
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
        precip = Precipitation1M(
            CMP.CloudLiquid(FT, toml_dict),
            CMP.Rain(FT, toml_dict),
            CMP.CollisionEff(FT, toml_dict),
            rf,
            CMP.Blk1MVelType(FT, toml_dict),
        )
    elseif precipitation_choice == "Precipitation2M"
        if rain_formation_choice == "SB2006"
            precip = Precipitation2M(CMP.SB2006(FT, toml_dict))
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
    else
        error("Invalid precipitation choice: $precipitation_choice")
    end
    return precip
end
