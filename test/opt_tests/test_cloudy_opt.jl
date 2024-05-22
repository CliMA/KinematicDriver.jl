

include(joinpath(pkgdir(KinematicDriver), "test", "create_parameters.jl"))
include(joinpath(pkgdir(KinematicDriver), "test", "plotting_utils.jl"))

function test_cloudy_allocation(::Type{FT}) where {FT}

    # setup
    default_toml_dict = CP.create_toml_dict(FT)
    toml_dict = override_toml_dict("_", default_toml_dict)
    common_params = create_common_parameters(toml_dict)
    kid_params = create_kid_parameters(toml_dict)
    thermo_params = create_thermodynamics_parameters(toml_dict)
    air_params = CMP.AirProperties(toml_dict)
    activation_params = CMP.AerosolActivationParameters(toml_dict)

    moisture = CO.get_moisture_type("CloudyMoisture", toml_dict)
    precip = CO.get_precipitation_type("CloudyPrecip", toml_dict)
    TS = CO.TimeStepping(FT(1), FT(10), FT(100))

    space, face_space = K1D.make_function_space(FT, FT(0), FT(3000), 64)
    coord = CC.Fields.coordinate_field(space)

    ρ_profile = CO.ρ_ivp(FT, kid_params, thermo_params)
    cloudy_params, cloudy_pdists = create_cloudy_parameters(FT)
    init = map(
        coord -> CO.cloudy_initial_condition(
            cloudy_pdists,
            CO.initial_condition_1d(FT, common_params, kid_params, thermo_params, ρ_profile, coord.z),
        ),
        coord,
    )

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

    # test
    @test_opt CO.precompute_aux_thermo!(moisture, Y, aux)
    CO.precompute_aux_thermo!(moisture, Y, aux)
    @test 64 >= @allocated CO.precompute_aux_thermo!(moisture, Y, aux)

    @test_opt CO.precompute_aux_precip!(precip, Y, aux)
    CO.precompute_aux_precip!(precip, Y, aux)
    @test 944 >= @allocated CO.precompute_aux_precip!(precip, Y, aux)

    @test_opt CO.precompute_aux_moisture_sources!(moisture, dY, Y, aux, 0.0)
    CO.precompute_aux_moisture_sources!(moisture, dY, Y, aux, 0.0)
    @test 0 == @allocated CO.precompute_aux_moisture_sources!(moisture, dY, Y, aux, 0.0)

    @test_opt CO.precompute_aux_precip_sources!(precip, aux)
    CO.precompute_aux_precip_sources!(precip, aux)
    @test 928 >= @allocated CO.precompute_aux_precip_sources!(precip, aux)

    @test_opt K1D.precompute_aux_activation!(precip, dY, Y, aux, 0.0)
    K1D.precompute_aux_activation!(precip, dY, Y, aux, 0.0)
    @test 624 >= @allocated K1D.precompute_aux_activation!(precip, dY, Y, aux, 0.0)

end

test_cloudy_allocation(Float64)
