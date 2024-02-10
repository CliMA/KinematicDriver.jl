"""
    Create Box model parameters.
"""
module Parameters

abstract type AbstractKinematicParameters end
const AKP = AbstractKinematicParameters

# Define KiD parameters
Base.@kwdef struct BoxModelParameters{FT} <: AKP
    ρ_air::FT
    precip_sources::Int
    precip_sinks::Int
    prescribed_Nd::FT
end

ρ_air(ps::AKP) = ps.ρ_air
precip_sources(ps::AKP) = ps.precip_sources
precip_sinks(ps::AKP) = ps.precip_sinks
prescribed_Nd(ps::AKP) = ps.prescribed_Nd

Base.eltype(::BoxModelParameters{FT}) where {FT} = FT
# Magic needed to get rid of length(ps) error
Base.broadcastable(ps::AKP) = Ref(ps)

end
