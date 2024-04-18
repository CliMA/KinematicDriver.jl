module Parameters

using DocStringExtensions
import ClimaParams as CP

abstract type AbstractKinematicParameters end
const AKP = AbstractKinematicParameters

"""
    Kinematic1DParameters{FT}

    Free parameters for the kinematic 1-dimensional simulation

    #Fields
    $(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Kinematic1DParameters{FT} <: AKP
    "Maximum updraft momentum flux m/s kg/m3"
    w1::FT
    "Time when the updraft is switched off [s]"
    t1::FT
    "Surface pressure [Pa]"
    p0::FT
    "Switch to include flux correction for moisture transport"
    qtot_flux_correction::Int
    "Assumed aerosol dry radius [m]"
    r_dry::FT
    "Assumed aerosol standard deviation [-]"
    std_dry::FT
    "Assumed aerosol hygroscopicity"
    κ::FT
end

function Kinematic1DParameters(td::CP.AbstractTOMLDict)
    #name_map = (;
    #    :prescribed_flow_w1 => :w1,
    #    :prescribed_flow_t1 => :t1,
    #    :surface_pressure => :p0,
    #    :qtot_flux_correction_flag => :qtot_flux_correction,
    #    :r_dry => :r_dry,
    #    :std_dry => :std_dry,
    #    :kappa => :κ,
    #)
    #parameters = CP.get_parameter_values(td, name_map, "Kinematic1D")
    #FT = CP.float_type(td)
    #return Kinematic1DParameters{FT}(; parameters)
    (; data) = td
    FT = CP.float_type(td)
    return Kinematic1DParameters(
        FT(data["prescribed_flow_w1"]["value"]),
        FT(data["prescribed_flow_t1"]["value"]),
        FT(data["surface_pressure"]["value"]),
        Int(data["qtot_flux_correction_flag"]["value"]),
        FT(data["r_dry"]["value"]),
        FT(data["std_dry"]["value"]),
        FT(data["kappa"]["value"]),
    )
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
