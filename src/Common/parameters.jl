"""
    Create common model parameters.

    (Overwriting the default values happens in the driver file.
    Here we are just creating the CommonParameters struct)
"""
module Parameters

abstract type AbstractCommonParameters end
const ACP = AbstractCommonParameters

# Define KiD parameters
Base.@kwdef struct CommonParameters{FT} <: ACP
    precip_sources::Int
    precip_sinks::Int
    prescribed_Nd::FT
end

precip_sources(ps::ACP) = ps.precip_sources
precip_sinks(ps::ACP) = ps.precip_sinks
prescribed_Nd(ps::ACP) = ps.prescribed_Nd

Base.eltype(::CommonParameters{FT}) where {FT} = FT
# Magic needed to get rid of length(ps) error
Base.broadcastable(ps::ACP) = Ref(ps)
end
