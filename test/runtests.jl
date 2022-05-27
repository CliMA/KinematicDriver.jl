using Test

import Interpolations

import CLIMAParameters
import ClimaCore

include("../src/Kinematic1D.jl")

const IP = Interpolations
const CP = CLIMAParameters
const CC = ClimaCore

const KiD = Kinematic1D

# Instantiate CliMA Parameters and overwrite the defaults to match PySDM
params = KiD.params_overwrite

include("unit_test.jl")
include("initial_profile_test.jl")
