module Parameters

using DocStringExtensions
import ClimaParams as CP
import Cloudy as CL

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
    "Switch to define aerosol activation approach; open system vs. closed system"
    open_system_activation::Bool
    "Switch to define aerosol activation approach; local vs. cloud-base"
    local_activation::Bool
end

function CommonParameters(td::CP.AbstractTOMLDict)
    name_map = (;
        :precipitation_sources_flag => :precip_sources,
        :precipitation_sinks_flag => :precip_sinks,
        :prescribed_Nd => :prescribed_Nd,
        :open_system_activation => :open_system_activation,
        :local_activation => :local_activation,
    )
    parameters = CP.get_parameter_values(td, name_map, "Common")
    FT = CP.float_type(td)
    return CommonParameters{FT}(; parameters...)
end

"""
    CloudyParameters{FT}

    Parameters for the kinematic driver simulations using Cloudy

    #Fields
    $(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct CloudyParameters{FT, NM, ND, N, P, T} <: ACP
    "Number of prognostic moments associated with each dist"
    NProgMoms::NTuple{ND, Int}
    "Normalizing number density and mass"
    norms::Tuple{FT, FT}
    "Normalizing moments"
    mom_norms::NTuple{NM, FT}
    "Coalescence Data"
    coal_data::CL.Coalescence.CoalescenceData{N, P, FT, T}
    "Vel information for power law"
    vel::Tuple{Tuple{FT, FT}}
    "Size threshold for distinguishing cloud and rain"
    size_threshold::FT

    function CloudyParameters(
        NProgMoms::NTuple{ND, Int},
        norms::Tuple{FT, FT},
        mom_norms::NTuple{NM, FT},
        coal_data::CL.Coalescence.CoalescenceData{N, P, FT, T},
        vel::Tuple{Tuple{FT, FT}},
        size_threshold::FT = 5e-10,
    ) where {NM, ND, N, P, T, FT <: Real}
        normed_vel = map(vel) do v
            (v[1] * norms[2]^v[2], v[2])
        end
        new{FT, NM, ND, N, P, T}(NProgMoms, norms, mom_norms, coal_data, normed_vel, size_threshold)
    end
end

precip_sources(ps::ACP) = ps.precip_sources
precip_sinks(ps::ACP) = ps.precip_sinks
prescribed_Nd(ps::ACP) = ps.prescribed_Nd

Base.eltype(::CommonParameters{FT}) where {FT} = FT
Base.eltype(::CloudyParameters{FT}) where {FT} = FT
# Magic needed to get rid of length(ps) error
Base.broadcastable(ps::ACP) = Ref(ps)
Base.broadcastable(x::CommonParameters) = Ref(x)
Base.broadcastable(x::CloudyParameters) = Ref(x)

end
