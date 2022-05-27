using Test

import Interpolations
const IP = Interpolations

import CLIMAParameters
const CP = CLIMAParameters

import ClimaCore
const CC = ClimaCore

include("../src/Kinematic1D.jl")
const KiD = Kinematic1D

# Instantiate CliMA Parameters and overwrite the defaults to match PySDM
# TODO - move it to InternalParameters (similar as in TC or CM)
struct EarthParameterSet{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end
CP.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
CP.gas_constant(ps::EarthParameterSet) = ps.nt.gas_constant
CP.Planet.molmass_dryair(ps::EarthParameterSet) = ps.nt.molmass_dryair
CP.Planet.cp_d(ps::EarthParameterSet) = ps.nt.cp_d
nt = (; MSLP = 100000.0, gas_constant = 8.314462618, molmass_dryair = 0.02896998, cp_d = 1005.0)
params = EarthParameterSet(nt)

include("unit_test.jl")
include("initial_profile_test.jl")
