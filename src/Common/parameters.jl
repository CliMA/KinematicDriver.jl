module Parameters

using DocStringExtensions
import ClimaParams as CP

abstract type AbstractCommonParameters end
const ACP = AbstractCommonParameters

"""
    CommonParameters{FT}

    Common parameters for the kinematic driver simulations

    #Fields
    $(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct CommonParameters{FT} <: ACP
    "Switch to include precipitation formation terms"
    precip_sources::Bool
    "Switch to include precipitation sink terms"
    precip_sinks::Bool
    "Prescribed number concentration of cloud droplets (needed for some schemes) [1/m3]"
    prescribed_Nd::FT
end

function CommonParameters(td::CP.AbstractTOMLDict)
    name_map = (;
        :precipitation_sources_flag => :precip_sources,
        :precipitation_sinks_flag => :precip_sinks,
        :prescribed_Nd => :prescribed_Nd,
    )
    parameters = CP.get_parameter_values(td, name_map, "Common")
    FT = CP.float_type(td)
    return CommonParameters{FT}(; parameters...)
end

precip_sources(ps::ACP) = ps.precip_sources
precip_sinks(ps::ACP) = ps.precip_sinks
prescribed_Nd(ps::ACP) = ps.prescribed_Nd

Base.eltype(::CommonParameters{FT}) where {FT} = FT
# Magic needed to get rid of length(ps) error
Base.broadcastable(ps::ACP) = Ref(ps)
Base.broadcastable(x::CommonParameters) = Ref(x)

end
