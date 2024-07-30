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
export MoistureP3
export CloudyMoisture

export NoPrecipitation
export Precipitation0M
export Precipitation1M
export Precipitation2M
export PrecipitationP3
export CloudyPrecip

abstract type AbstractStyle end
Base.broadcastable(x::AbstractStyle) = Ref(x)

abstract type AbstractMoistureStyle <: AbstractStyle end
abstract type AbstractPrecipitationStyle <: AbstractStyle end

#struct NoMoisture <: AbstractMoistureStyle end
struct CloudyMoisture <: AbstractMoistureStyle end
struct EquilibriumMoisture <: AbstractMoistureStyle end
struct NonEquilibriumMoisture{CL, IC} <: AbstractMoistureStyle
    liquid::CL
    ice::IC
end
struct MoistureP3 <: AbstractMoistureStyle end

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
struct PrecipitationP3{P3, CH, PT} <: AbstractPrecipitationStyle
    p3_params::P3
    Chen2022::CH
    sb2006::PT
end
struct CloudyPrecip <: AbstractPrecipitationStyle end
