"""
Types for dispatching on different equation sets for KiD
"""

export AbstractStyle
export AbstractMoistureStyle
export AbstractPrecipitationStyle

# TODO - add support for a dry case
#export NoMoisture
export EquilibriumMoisture
export NonEquilibriumMoisture
export EquilibriumMoisture_ρθq
export EquilibriumMoisture_ρdTq
export NonEquilibriumMoisture_ρθq
export NonEquilibriumMoisture_ρdTq

export NoPrecipitation
export Precipitation0M
export Precipitation1M
export Precipitation2M

abstract type AbstractStyle end
Base.broadcastable(x::AbstractStyle) = Ref(x)

abstract type AbstractMoistureStyle <: AbstractStyle end
abstract type AbstractPrecipitationStyle <: AbstractStyle end

#struct NoMoisture <: AbstractMoistureStyle end
abstract type EquilibriumMoisture <: AbstractMoistureStyle end
abstract type NonEquilibriumMoisture <: AbstractMoistureStyle end
struct EquilibriumMoisture_ρθq <: EquilibriumMoisture end
struct EquilibriumMoisture_ρdTq <: EquilibriumMoisture end
struct NonEquilibriumMoisture_ρθq{CL, IC} <: NonEquilibriumMoisture
    liquid::CL
    ice::IC
end
struct NonEquilibriumMoisture_ρdTq{CL, IC} <: NonEquilibriumMoisture
    liquid::CL
    ice::IC
end

struct NoPrecipitation <: AbstractPrecipitationStyle end
struct Precipitation0M{P0M} <: AbstractPrecipitationStyle
    params::P0M
end
struct Precipitation1M{CL, IC, RN, SN, CE, PT, ST} <: AbstractPrecipitationStyle
    liquid::CL
    ice::IC
    rain::RN
    snow::SN
    ce::CE
    rain_formation::PT
    sedimentation::ST
end
struct Precipitation2M{PT, ST} <: AbstractPrecipitationStyle
    rain_formation::PT
    sedimentation::ST
end
