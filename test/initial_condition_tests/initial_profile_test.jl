"""
Test that initial profiles match between CliMA and PySDM
"""

using Test

import Interpolations as IP
import LinearAlgebra

import ClimaParams
import ClimaCore as CC
import Thermodynamics
import CloudMicrophysics.Parameters as CMP
import Kinematic1D
import Kinematic1D.K1DModel as KID

const kid_dir = pkgdir(Kinematic1D)
include(joinpath(kid_dir, "test", "data_utils.jl"))
include(joinpath(kid_dir, "test", "plotting_utils.jl"))
include(joinpath(kid_dir, "test", "create_parameters.jl"))

const FT = Float64

# Create parameters overwrite the defaults to match PySDM
default_toml_dict = CP.create_toml_dict(FT)
toml_dict = override_toml_dict(@__DIR__, default_toml_dict)
thermo_params = create_thermodynamics_parameters(toml_dict)
common_params = create_common_parameters(toml_dict)
kid_params = create_kid_parameters(toml_dict)
air_params = CMP.AirProperties(FT, toml_dict)
activation_params = CMP.AerosolActivationParameters(FT, toml_dict)
params = (common_params, kid_params, thermo_params, air_params, activation_params)

function compare_profiles(; is_dry_flag::Bool)
    # Computational domain ...
    z_min = FT(0)
    z_max = FT(3220)
    n_elem = 222
    # ... and the created coordinates
    space, face_space = KID.make_function_space(FT, z_min, z_max, n_elem)
    coord = CC.Fields.coordinate_field(space)
    face_coord = CC.Fields.coordinate_field(face_space)

    # Solve the initial value problem for density profile
    ρ_profile = CO.ρ_ivp(FT, kid_params, thermo_params, dry = is_dry_flag)
    # Create the initial condition profiles
    init = map(
        coord -> CO.initial_condition_1d(
            FT,
            common_params,
            kid_params,
            thermo_params,
            ρ_profile,
            coord.z,
            dry = is_dry_flag,
        ),
        coord,
    )

    # Store the CliMA KID initial profiles
    z_centers = parent(CC.Fields.coordinate_field(space))
    T = parent(init.T)
    p = parent(init.p)
    ρ = parent(init.ρ)
    θ_dry = parent(init.θ_dry)
    q_vap = parent(init.q_tot) - parent(init.q_liq) - parent(init.q_ice)
    q_liq = parent(init.q_liq)
    KM_data = (; z_centers, q_vap, ρ, θ_dry, T, p, q_liq)

    # Read in the PySDM KID initial profiles
    sdm_case = is_dry_flag ? "dry" : "wet"
    sdm_data = load_sdm_data(sdm_case)

    # Interpolate both (we are not doing computations on the same grid)
    KM_T = IP.linear_interpolation(KM_data.z_centers[:], KM_data.T[:], extrapolation_bc = IP.Line())
    KM_p = IP.linear_interpolation(KM_data.z_centers[:], KM_data.p[:], extrapolation_bc = IP.Line())
    KM_ρ = IP.linear_interpolation(KM_data.z_centers[:], KM_data.ρ[:], extrapolation_bc = IP.Line())
    KM_q_liq = IP.linear_interpolation(KM_data.z_centers[:], KM_data.q_liq[:], extrapolation_bc = IP.Line())
    KM_q_vap = IP.linear_interpolation(KM_data.z_centers[:], KM_data.q_vap[:], extrapolation_bc = IP.Line())
    KM_θ_dry = IP.linear_interpolation(KM_data.z_centers[:], KM_data.θ_dry[:], extrapolation_bc = IP.Line())

    SD_T = IP.linear_interpolation(sdm_data.z_sdm, sdm_data.T_sdm, extrapolation_bc = IP.Line())
    SD_p = IP.linear_interpolation(sdm_data.z_sdm, sdm_data.P_sdm, extrapolation_bc = IP.Line())
    SD_ρ = IP.linear_interpolation(sdm_data.z_sdm, sdm_data.rho_sdm, extrapolation_bc = IP.Line())
    SD_q_liq = IP.linear_interpolation(sdm_data.z_sdm, sdm_data.ql_sdm, extrapolation_bc = IP.Line())
    SD_q_vap = IP.linear_interpolation(sdm_data.z_sdm, sdm_data.qv_sdm, extrapolation_bc = IP.Line())
    SD_θ_dry = IP.linear_interpolation(sdm_data.z_sdm, sdm_data.thetad_sdm, extrapolation_bc = IP.Line())

    # Test that the initial profiles match between CliMA and PySDM
    test_title = "Case: " * sdm_case * ". Initial profiles of water (vapour and liquid):"
    if is_dry_flag
        @testset "$test_title" begin
            z_test = z_min:100:z_max
            @test all(isapprox(KM_q_vap(z_test), SD_q_vap(z_test), rtol = 1e-6))
            @test all(isapprox(KM_q_liq(z_test), SD_q_liq(z_test), rtol = 1e-6))
        end
    else
        @testset "$test_title" begin
            z_test = z_min:100:z_max
            @test all(isapprox(KM_q_vap(z_test), SD_q_vap(z_test), rtol = 1e-4))
            @test all(isapprox(KM_q_liq(z_test), SD_q_liq(z_test), rtol = 1e-6))
        end
    end
    test_title = "Case: " * sdm_case * ". Initial profiles of (T, p, ρ):"
    @testset "$test_title" begin
        z_test = z_min:100:z_max
        @test all(isapprox(KM_T(z_test), SD_T(z_test), rtol = 1e-2))
        @test all(isapprox(KM_p(z_test), SD_p(z_test), rtol = 1e-2))
        @test all(isapprox(KM_ρ(z_test), SD_ρ(z_test), rtol = 1e-2))
    end

    test_title = "Case: " * sdm_case * ". Initial profiles of θ_dry:"
    @testset "$test_title" begin
        z_test = z_min:100:z_max
        @test all(isapprox(KM_θ_dry(z_test), SD_θ_dry(z_test), rtol = 1e-6))
    end

    plot_initial_profiles_comparison(KM_data, sdm_case = sdm_case)
end

compare_profiles(is_dry_flag = true)
compare_profiles(is_dry_flag = false)
