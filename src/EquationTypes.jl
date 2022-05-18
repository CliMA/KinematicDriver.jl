"""
Types for dispatching on different equation sets for KiD
"""

export AbstractStyle
export AbstractMoistureStyle
export AbstractPrecipitationStyle

export NoMoisture
export EquilibriumMoisture
export NonEquilibriumMoisture

export NoPrecipitation
export Precipitation0M
export Precipitation1M

abstract type AbstractStyle end
abstract type AbstractMoistureStyle <: AbstractStyle end
abstract type AbstractPrecipitationStyle <: AbstractStyle end

struct NoMoisture <: AbstractMoistureStyle end
struct EquilibriumMoisture <: AbstractMoistureStyle end
struct NonEquilibriumMoisture <: AbstractMoistureStyle end

struct NoPrecipitation <: AbstractPrecipitationStyle end
struct Precipitation0M <: AbstractPrecipitationStyle end
struct Precipitation1M <: AbstractPrecipitationStyle end
