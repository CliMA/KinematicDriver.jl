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

include(joinpath(pkgdir(KinematicDriver), "test", "create_parameters.jl"))
include(joinpath(pkgdir(KinematicDriver), "test", "plotting_utils.jl"))

function get_tendency_function_arguments(::Type{FT}, moisture_choice, precipitation_choice) where {FT}

    # setup
    default_toml_dict = CP.create_toml_dict(FT)
    toml_dict = override_toml_dict("_", default_toml_dict)
    common_params = create_common_parameters(toml_dict)
    kid_params = create_kid_parameters(toml_dict)
    thermo_params = create_thermodynamics_parameters(toml_dict)
    air_params = CMP.AirProperties(toml_dict)
    activation_params = CMP.AerosolActivationParameters(toml_dict)

    moisture = CO.get_moisture_type(moisture_choice, toml_dict)
    precip = CO.get_precipitation_type(precipitation_choice, toml_dict)
    TS = CO.TimeStepping(FT(1), FT(10), FT(100))

    space, face_space = K1D.make_function_space(FT, FT(0), FT(3000), 64)
    coord = CC.Fields.coordinate_field(space)

    ρ_profile = CO.ρ_ivp(FT, kid_params, thermo_params)
    if precipitation_choice == "CloudyPrecip"
        cloudy_params, cloudy_pdists = create_cloudy_parameters(FT)
        init = map(
            coord -> CO.cloudy_initial_condition(
                cloudy_pdists,
                CO.initial_condition_1d(FT, common_params, kid_params, thermo_params, ρ_profile, coord.z),
            ),
            coord,
        )
    else
        cloudy_params = nothing
        init = map(
            coord -> CO.initial_condition_1d(FT, common_params, kid_params, thermo_params, ρ_profile, coord.z),
            coord,
        )
    end

    aux = K1D.initialise_aux(
        FT,
        init,
        common_params,
        kid_params,
        thermo_params,
        air_params,
        activation_params,
        TS,
        nothing,
        face_space,
        moisture,
        precip,
        cloudy_params,
    )

    Y = CO.initialise_state(moisture, precip, init)
    dY = Y

    return moisture, precip, Y, dY, aux
end

@testset "NonEquilibriumMoisture + Precipitation2M optimization tests" begin
    # setup
    moisture, precip, Y, dY, aux = get_tendency_function_arguments(Float64, "NonEquilibriumMoisture", "Precipitation2M")

    # test
    @test_opt CO.precompute_aux_thermo!(moisture, Y, aux)
    CO.precompute_aux_thermo!(moisture, Y, aux)
    @test 64 >= @allocated CO.precompute_aux_thermo!(moisture, Y, aux)

    @test_opt CO.precompute_aux_precip!(precip, Y, aux)
    CO.precompute_aux_precip!(precip, Y, aux)
    @test 32 >= @allocated CO.precompute_aux_precip!(precip, Y, aux)

    # @test_opt CO.precompute_aux_moisture_sources!(moisture, aux) #TODO
    CO.precompute_aux_moisture_sources!(moisture, aux)
    @test 7200 == @allocated CO.precompute_aux_moisture_sources!(moisture, aux)

    # @test_opt CO.precompute_aux_precip_sources!(precip, aux) #TODO
    CO.precompute_aux_precip_sources!(precip, aux)
    @test 153000 >= @allocated CO.precompute_aux_precip_sources!(precip, aux)

    @test_opt K1D.precompute_aux_activation!(precip, dY, Y, aux, 0.0)
    K1D.precompute_aux_activation!(precip, dY, Y, aux, 0.0)
    # TODO allocation occurs for finding cloud base
    @test 9824 >= @allocated K1D.precompute_aux_activation!(precip, dY, Y, aux, 0.0)
end

@testset "Cloudy optimization tests" begin
    # setup
    moisture, precip, Y, dY, aux = get_tendency_function_arguments(Float64, "CloudyMoisture", "CloudyPrecip")

    # test
    @test_opt CO.precompute_aux_thermo!(moisture, Y, aux)
    CO.precompute_aux_thermo!(moisture, Y, aux)
    @test 64 >= @allocated CO.precompute_aux_thermo!(moisture, Y, aux)

    @test_opt CO.precompute_aux_precip!(precip, Y, aux)
    CO.precompute_aux_precip!(precip, Y, aux)
    @test 928 >= @allocated CO.precompute_aux_precip!(precip, Y, aux)

    @test_opt CO.precompute_aux_moisture_sources!(moisture, aux)
    CO.precompute_aux_moisture_sources!(moisture, aux)
    @test 0 == @allocated CO.precompute_aux_moisture_sources!(moisture, aux)

    @test_opt CO.precompute_aux_precip_sources!(precip, aux)
    CO.precompute_aux_precip_sources!(precip, aux)
    @test 928 >= @allocated CO.precompute_aux_precip_sources!(precip, aux)

    @test_opt K1D.precompute_aux_activation!(precip, dY, Y, aux, 0.0)
    K1D.precompute_aux_activation!(precip, dY, Y, aux, 0.0)
    # TODO allocation occurs for finding cloud base
    @test 10240 >= @allocated K1D.precompute_aux_activation!(precip, dY, Y, aux, 0.0)
end
