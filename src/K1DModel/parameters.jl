"""
    Create Kinematic1D model parameters.

    (Overwriting the default values happens in the driver file.
    Here we are just creating the Kinematic1DParameters struct)
"""
module Parameters

abstract type AbstractKinematicParameters end
const AKP = AbstractKinematicParameters

# Define KiD parameters
Base.@kwdef struct Kinematic1DParameters{FT} <: AKP
    w1::FT
    t1::FT
    p0::FT
    qtot_flux_correction::Int
    r_dry::FT
    std_dry::FT
    κ::FT
end

w1(ps::AKP) = ps.w1
t1(ps::AKP) = ps.t1
p0(ps::AKP) = ps.p0
#aerosol parameters for 2M scheme
r_dry(ps::AKP) = ps.r_dry
std_dry(ps::AKP) = ps.std_dry
κ(ps::AKP) = ps.κ

Base.eltype(::Kinematic1DParameters{FT}) where {FT} = FT
# Magic needed to get rid of length(ps) error
Base.broadcastable(ps::AKP) = Ref(ps)
end
