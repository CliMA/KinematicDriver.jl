"""
    Create Kinematic1D model parameters.

    (Overwriting the default values happens in the driver file.
    Here we are just creating the KinematicParameters struct)
"""
module Parameters

abstract type AbstractKinematicParameters end
const AKP = AbstractKinematicParameters

# Define KiD parameters
Base.@kwdef struct KinematicParameters{FT} <: AKP
    w1::FT
    t1::FT
    p0::FT
    precip_sources::Int
    precip_sinks::Int
    prescribed_Nd::FT
    qtot_flux_correction::Int
    r_dry::FT
    std_dry::FT
    κ::FT
end

w1(ps::AKP) = ps.w1
t1(ps::AKP) = ps.t1
precip_sources(ps::AKP) = ps.precip_sources
precip_sinks(ps::AKP) = ps.precip_sinks
prescribed_Nd(ps::AKP) = ps.prescribed_Nd
#aerosol parameters for 2M scheme
r_dry(ps::AKP) = ps.r_dry
std_dry(ps::AKP) = ps.std_dry
κ(ps::AKP) = ps.κ

Base.eltype(::KinematicParameters{FT}) where {FT} = FT
# Magic needed to get rid of length(ps) error
Base.broadcastable(ps::AKP) = Ref(ps)
end
