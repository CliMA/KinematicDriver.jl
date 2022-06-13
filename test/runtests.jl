using Test

import Interpolations
import LinearAlgebra

import CLIMAParameters
import ClimaCore
import Thermodynamics

include("../src/Kinematic1D.jl")

const IP = Interpolations
const LA = LinearAlgebra
const CP = CLIMAParameters
const CC = ClimaCore
const TD = Thermodynamics

const KiD = Kinematic1D

# Instantiate CliMA Parameters and overwrite the defaults to match PySDM
params = KiD.params_overwrite

include("unit_test.jl")
include("initial_profile_test.jl")
