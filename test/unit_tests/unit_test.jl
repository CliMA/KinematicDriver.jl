"""
    Some elementary unit tests
"""

using Test
import LinearAlgebra as LA

import ClimaParams as CP
import ClimaCore as CC
import Thermodynamics as TD
import KinematicDriver
import KinematicDriver.Common as CO
import KinematicDriver.BoxModel as BX
import KinematicDriver.K1DModel as K1D
import KinematicDriver.K2DModel as K2D
import CloudMicrophysics.Parameters as CMP

const FT = Float64

const kid_dir = pkgdir(KinematicDriver)
include(joinpath(kid_dir, "test", "create_parameters.jl"))

# override the defaults
default_toml_dict = CP.create_toml_dict(FT)
toml_dict = override_toml_dict(@__DIR__, default_toml_dict)

# create all the parameters structs ...
common_params = create_common_parameters(toml_dict)
thermo_params = create_thermodynamics_parameters(toml_dict)
kid_params = create_kid_parameters(toml_dict)
air_params = CMP.AirProperties(toml_dict)
activation_params = CMP.AerosolActivationParameters(toml_dict)
# ... for cloud condensate options ...
equil_moist = CO.EquilibriumMoisture()
nequil_moist = CO.NonEquilibriumMoisture(CMP.CloudLiquid(toml_dict), CMP.CloudIce(toml_dict))
# ... and precipitation options
no_precip = CO.NoPrecipitation()
precip_0m = CO.Precipitation0M(CMP.Parameters0M(toml_dict))
precip_1m = CO.Precipitation1M(
    CMP.CloudLiquid(toml_dict),
    CMP.CloudIce(toml_dict),
    CMP.Rain(toml_dict),
    CMP.Snow(toml_dict),
    CMP.CollisionEff(toml_dict),
    CMP.KK2000(toml_dict),
    CMP.Blk1MVelType(toml_dict),
)
precip_2m = CO.Precipitation2M(CMP.SB2006(toml_dict), CMP.SB2006VelType(toml_dict))

# common unit tests
include("./common_unit_test.jl")

# KinematricDriver 1D unit tests
include("./k1d_unit_test.jl")

# KinematricDriver 2D unit tests
include("./k2d_unit_test.jl")

# calibration pipeline unit tests
include("./calibration_pipeline/run_calibration_pipeline_unit_tests.jl")
