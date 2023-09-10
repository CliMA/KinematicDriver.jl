"""
    Interface with free parameters needed by KiD, CloudMicrophysics and Thermodynamics

    Overwriting the default values happens in the driver file.
    Here we are just creating the KinematicParameters struct and
    forwarding the parameters definitions to Thermodynamics and CloudMicrophysics definitions.
"""
module Parameters

import Thermodynamics as TD
import CloudMicrophysics as CM

abstract type AbstractKinematicParameters end
const AKP = AbstractKinematicParameters

# Define KiD parameters
Base.@kwdef struct KinematicParameters{FT, MP} <: AKP
    w1::FT
    t1::FT
    p0::FT
    z_0::FT
    z_1::FT
    z_2::FT
    rv_0::FT
    rv_1::FT
    rv_2::FT
    tht_0::FT
    tht_1::FT
    tht_2::FT
    precip_sources::Int
    precip_sinks::Int
    prescribed_Nd::FT
    qtot_flux_correction::Int
    r_dry::FT
    std_dry::FT
    κ::FT
    microphys_params::MP
end
thermodynamics_params(ps::AKP) = CM.Parameters.thermodynamics_params(ps.microphys_params)
microphysics_params(ps::AKP) = ps.microphys_params

Base.eltype(::KinematicParameters{FT}) where {FT} = FT
w1(ps::AKP) = ps.w1
t1(ps::AKP) = ps.t1

z_0(ps::AKP) = ps.z_0
z_1(ps::AKP) = ps.z_1
z_2(ps::AKP) = ps.z_2
rv_0(ps::AKP) = ps.rv_0
rv_1(ps::AKP) = ps.rv_1
rv_2(ps::AKP) = ps.rv_2
tht_0(ps::AKP) = ps.tht_0
tht_1(ps::AKP) = ps.tht_1
tht_2(ps::AKP) = ps.tht_2
precip_sources(ps::AKP) = ps.precip_sources
precip_sinks(ps::AKP) = ps.precip_sinks
prescribed_Nd(ps::AKP) = ps.prescribed_Nd
#aerosol parameters for 2M scheme
r_dry(ps::AKP) = ps.r_dry
std_dry(ps::AKP) = ps.std_dry
κ(ps::AKP) = ps.κ

# Forward parameters to Thermodynamics
const TDPS = TD.Parameters.ThermodynamicsParameters
for var in fieldnames(TDPS)
    @eval $var(ps::AKP) = TD.Parameters.$var(thermodynamics_params(ps))
end

# Define derived parameters
molmass_ratio(ps::AKP) = TD.Parameters.molmass_ratio(thermodynamics_params(ps))
R_d(ps::AKP) = TD.Parameters.R_d(thermodynamics_params(ps))
R_v(ps::AKP) = TD.Parameters.R_v(thermodynamics_params(ps))
cp_d(ps::AKP) = TD.Parameters.cp_d(thermodynamics_params(ps))
cv_l(ps::AKP) = TD.Parameters.cv_l(thermodynamics_params(ps))

# Magic needed to get rid of length(ps) error
Base.broadcastable(ps::AKP) = Ref(ps)
end
