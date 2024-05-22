"""
    Optimization tests
"""

import ClimaCore as CC
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics as CM
import KinematicDriver
import KinematicDriver.Common as CO
import KinematicDriver.K1DModel as K1D
using JET
using Test

# common unit tests
include("./test_cloudy_opt.jl")
