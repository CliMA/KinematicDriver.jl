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

const kid_dir = pkgdir(Kinematic1D)
include(joinpath(kid_dir, "test", "create_parameters.jl"))

# override the defaults
default_toml_dict = CP.create_toml_dict(FT, dict_type = "alias")
toml_dict = override_toml_dict(@__DIR__, default_toml_dict)

# create all the parameters structs ...
common_params = create_common_parameters(toml_dict)
thermo_params = create_thermodynamics_parameters(toml_dict)
kid_params = create_kid_parameters(toml_dict)
air_params = CMP.AirProperties(FT, toml_dict)
activation_params = CMP.AerosolActivationParameters(FT, toml_dict)
# ... for cloud condensate options ...
equil_moist_ρθq = CO.EquilibriumMoisture_ρθq()
equil_moist_ρdTq = CO.EquilibriumMoisture_ρdTq()
nequil_moist_ρθq = CO.NonEquilibriumMoisture_ρθq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
nequil_moist_ρdTq = CO.NonEquilibriumMoisture_ρdTq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
# ... and precipitation options
no_precip = CO.NoPrecipitation()
precip_0m = CO.Precipitation0M(CMP.Parameters0M(FT, toml_dict))
precip_1m = CO.Precipitation1M(
    CMP.CloudLiquid(FT, toml_dict),
    CMP.CloudIce(FT, toml_dict),
    CMP.Rain(FT, toml_dict),
    CMP.Snow(FT, toml_dict),
    CMP.CollisionEff(FT, toml_dict),
    CMP.KK2000(FT, toml_dict),
    CMP.Blk1MVelType(FT, toml_dict),
)
precip_2m = CO.Precipitation2M(CMP.SB2006(FT, toml_dict), CMP.SB2006VelType(FT, toml_dict))

# common unit tests
include("./common_unit_test.jl")

# box model unit tests
include("./box_unit_test.jl")

# kinematic1d unit tests
include("./k1d_unit_test.jl")

# kinematic2d unit tests
include("./k2d_unit_test.jl")

# calibration pipeline unit tests
include("./calibration_pipeline/run_calibration_pipeline_unit_tests.jl")
