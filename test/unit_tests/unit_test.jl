"""
    Some elementary unit tests
"""

import Test as TT
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

FT = Float64

kid_dir = pkgdir(KinematicDriver)
include(joinpath(kid_dir, "test", "create_parameters.jl"))

# override the defaults
default_toml_dict = CP.create_toml_dict(FT)
toml_dict = override_toml_dict(default_toml_dict)

# create all the parameters structs ...
common_params = create_common_parameters(toml_dict)
thermo_params = create_thermodynamics_parameters(toml_dict)
kid_params = create_kid_parameters(toml_dict)
air_params = CMP.AirProperties(toml_dict)
activation_params = CMP.AerosolActivationParameters(toml_dict)
p3 = CMP.ParametersP3(FT)
Chen2022 = CMP.Chen2022VelType(FT)
# (use default boundary condition)
p3_boundary_condition = (;
    ice_start = false,
    _magnitude = Float64(0.5),
    _q_flux = Float64(0.65e-4),
    _N_flux = Float64(40000),
    _F_rim = Float64(0.2),
    _F_liq = Float64(0.2),
    _ρ_r_init = Float64(900),
)
# ... for cloud condensate options ...
equil_moist = CO.EquilibriumMoisture()
nequil_moist = CO.NonEquilibriumMoisture(CMP.CloudLiquid(toml_dict), CMP.CloudIce(toml_dict))
p3_moist = CO.MoistureP3()
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
precip_p3 = CO.PrecipitationP3(p3, Chen2022, CMP.SB2006(toml_dict), p3_boundary_condition)

TT.@testset "unit tests" begin
    # common unit tests
    include("./common_unit_test.jl")

    # KinematricDriver 1D unit tests
    include("./k1d_unit_test.jl")

    # KinematricDriver 2D unit tests
    include("./k2d_unit_test.jl")
end
nothing
