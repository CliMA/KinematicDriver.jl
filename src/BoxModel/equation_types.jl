"""
Types for dispatching on different equation sets
"""

export AbstractStyle
export AbstractPrecipitationStyle

export Precipitation1M
export Precipitation2M

abstract type AbstractStyle end
Base.broadcastable(x::AbstractStyle) = Ref(x)

abstract type AbstractPrecipitationStyle <: AbstractStyle end

struct Precipitation1M{CL, RN, CE, PT, ST} <: AbstractPrecipitationStyle
    liquid::CL
    rain::RN
    ce::CE
    rain_formation::PT
    sedimentation::ST
end
struct Precipitation2M{PT} <: AbstractPrecipitationStyle
    rain_formation::PT
end
