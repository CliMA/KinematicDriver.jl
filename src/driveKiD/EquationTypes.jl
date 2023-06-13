"""
Types for dispatching on different equation sets for KiD
"""

export AbstractStyle
export AbstractMoistureStyle
export AbstractPrecipitationStyle
export AbstractRainFormationStyle

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

export OneMomentRainFormation

abstract type AbstractStyle end
abstract type AbstractMoistureStyle <: AbstractStyle end
abstract type AbstractPrecipitationStyle <: AbstractStyle end
abstract type AbstractRainFormationStyle <: AbstractStyle end

Base.broadcastable(x::AbstractPrecipitationStyle) = Ref(x)

#struct NoMoisture <: AbstractMoistureStyle end
abstract type EquilibriumMoisture <: AbstractMoistureStyle end
abstract type NonEquilibriumMoisture <: AbstractMoistureStyle end
struct EquilibriumMoisture_ρθq <: EquilibriumMoisture end
struct EquilibriumMoisture_ρdTq <: EquilibriumMoisture end
struct NonEquilibriumMoisture_ρθq <: NonEquilibriumMoisture end
struct NonEquilibriumMoisture_ρdTq <: NonEquilibriumMoisture end

struct NoPrecipitation <: AbstractPrecipitationStyle end
struct Precipitation0M <: AbstractPrecipitationStyle end
struct Precipitation1M{PT} <: AbstractPrecipitationStyle
    rain_formation::PT
end
struct Precipitation2M{PT} <: AbstractPrecipitationStyle
    rain_formation::PT
end

struct OneMomentRainFormation <: AbstractRainFormationStyle end
