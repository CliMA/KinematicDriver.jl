"""
    Create Kinematic1D model parameters.

    (Overwriting the default values happens in the driver file.
    Here we are just creating the KinematicParameters struct)
"""
module Parameters

import CLIMAParameters as CP

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

function KinematicParameters(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return KinematicParameters(
        FT(data["prescribed_flow_w1"]["value"]),
        FT(data["prescribed_flow_t1"]["value"]),
        FT(data["surface_pressure"]["value"]),
        Int(data["precipitation_sources_flag"]["value"]),
        Int(data["precipitation_sinks_flag"]["value"]),
        FT(data["prescribed_Nd"]["value"]),
        Int(data["qtot_flux_correction_flag"]["value"]),
        FT(data["r_dry"]["value"]),
        FT(data["std_dry"]["value"]),
        FT(data["kappa"]["value"]),
    )
end
Base.eltype(::KinematicParameters{FT}) where {FT} = FT
# Magic needed to get rid of length(ps) error
Base.broadcastable(ps::AKP) = Ref(ps)

end #module
