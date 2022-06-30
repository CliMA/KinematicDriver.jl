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
    microphys_params::MP
end
thermodynamics_params(ps::AKP) = CM.Parameters.thermodynamics_params(ps.microphys_params)
microphysics_params(ps::AKP) = ps.microphys_params

Base.eltype(::KinematicParameters{FT}) where {FT} = FT
w1(ps::AKP) = ps.w1
t1(ps::AKP) = ps.t1

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
