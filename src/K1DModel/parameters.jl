module Parameters

using DocStringExtensions
import ClimaParams as CP
import ...Common as CO

abstract type AbstractKinematicParameters end
const AKP = AbstractKinematicParameters

"""
    ShipwayHill2012Parameters{FT}

    Free parameters for the kinematic 1-dimensional simulation

    #Fields
    $(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ShipwayHill2012Parameters{FT} <: CO.ShipwayHill2012
    "Maximum updraft momentum flux m/s kg/m3"
    w1::FT
    "Time when the updraft is switched off [s]"
    t1::FT
    "Surface pressure [Pa]"
    p0::FT
    "Initial profile parameter [m]"
    z_0::FT
    "Initial profile parameter [m]"
    z_1::FT
    "Initial profile parameter [m]"
    z_2::FT
    "Initial profile parameter [kg/kg]"
    rv_0::FT
    "Initial profile parameter [kg/kg]"
    rv_1::FT
    "Initial profile parameter [kg/kg]"
    rv_2::FT
    "Initial profile parameter [K]"
    tht_0::FT
    "Initial profile parameter [K]"
    tht_1::FT
    "Initial profile parameter [K]"
    tht_2::FT
    "Switch to include flux correction for moisture transport"
    qtot_flux_correction::Int
    "Assumed aerosol dry radius [m]"
    r_dry::FT
    "Assumed aerosol standard deviation [-]"
    std_dry::FT
    "Assumed aerosol hygroscopicity"
    κ::FT
end

function ShipwayHill2012Parameters(td::CP.AbstractTOMLDict)
    name_map = (;
        :prescribed_flow_w1 => :w1,
        :prescribed_flow_t1 => :t1,
        :surface_pressure => :p0,
        :init_cond_z0 => :z_0,
        :init_cond_z1 => :z_1,
        :init_cond_z2 => :z_2,
        :init_cond_rv0 => :rv_0,
        :init_cond_rv1 => :rv_1,
        :init_cond_rv2 => :rv_2,
        :init_cond_theta0 => :tht_0,
        :init_cond_theta1 => :tht_1,
        :init_cond_theta2 => :tht_2,
        :qtot_flux_correction_flag => :qtot_flux_correction,
        :r_dry => :r_dry,
        :std_dry => :std_dry,
        :kappa => :κ,
    )
    parameters = CP.get_parameter_values(td, name_map, "KinematicDriver")
    FT = CP.float_type(td)
    return ShipwayHill2012Parameters{FT}(; parameters...)
end

"""
    Jouan2020Parameters{FT}

    Free parameters for the kinematic 1-dimensional simulation

    #Fields
    $(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Jouan2020Parameters{FT} <: CO.Jouan2020
    "Maximum updraft momentum flux m/s kg/m3"
    w1::FT
    "Time when the updraft velocty switches sign [s]"
    t1::FT
    "Switch to include flux correction for moisture transport"
    qtot_flux_correction::Int
    "Assumed aerosol dry radius [m]"
    r_dry::FT
    "Assumed aerosol standard deviation [-]"
    std_dry::FT
    "Assumed aerosol hygroscopicity"
    κ::FT
end

function Jouan2020Parameters(td::CP.AbstractTOMLDict)
    name_map = (;
        :prescribed_flow_w1 => :w1,
        :prescribed_flow_t1 => :t1,
        :qtot_flux_correction_flag => :qtot_flux_correction,
        :r_dry => :r_dry,
        :std_dry => :std_dry,
        :kappa => :κ,
    )
    parameters = CP.get_parameter_values(td, name_map, "KinematicDriver")
    FT = CP.float_type(td)
    return Jouan2020Parameters{FT}(; parameters...)
end

# wrappers
for fn in fieldnames(ShipwayHill2012Parameters)
    @eval $(fn)(ps::ShipwayHill2012Parameters) = ps.$(fn)
end
for fn in fieldnames(Jouan2020Parameters)
    @eval $(fn)(ps::Jouan2020Parameters) = ps.$(fn)
end

Base.eltype(::ShipwayHill2012Parameters{FT}) where {FT} = FT
Base.eltype(::Jouan2020Parameters{FT}) where {FT} = FT

# Magic needed to get rid of length(ps) error
Base.broadcastable(ps::AKP) = Ref(ps)
Base.broadcastable(x::ShipwayHill2012Parameters) = Ref(x)
Base.broadcastable(x::Jouan2020Parameters) = Ref(x)

end
