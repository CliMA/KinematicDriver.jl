"""
    Some elementary unit tests
"""

using Test
import LinearAlgebra as LA

import CLIMAParameters as CP
import ClimaCore as CC
import Thermodynamics as TD
import Kinematic1D
import Kinematic1D.Common as CO
import Kinematic1D.BoxModel as BX
import Kinematic1D.K1DModel as K1D
import Kinematic1D.K2DModel as K2D
import CloudMicrophysics.Parameters as CMP

const FT = Float64

# box model unit tests
include("./box_unit_test.jl")

# kinematic1d unit tests
include("./k1d_unit_test.jl")

# kinematic2d unit tests
include("./k2d_unit_test.jl")

# calibration pipeline unit tests
include("./calibration_pipeline/run_calibration_pipeline_unit_tests.jl")
