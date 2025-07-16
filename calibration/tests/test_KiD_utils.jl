TT.@testset "creating parameters object" begin
    #setup
    u = [1e-4, 12345.0]
    u_names = ["cloud_liquid_water_specific_humidity_autoconversion_threshold", "rain_autoconversion_timescale"]
    model_settings = get_model_config()
    #action
    update_parameters!(model_settings, u, u_names)
    #test
    TT.@test model_settings["toml_dict"]["cloud_liquid_water_specific_humidity_autoconversion_threshold"]["value"] ==
             1e-4
    TT.@test model_settings["toml_dict"]["rain_autoconversion_timescale"]["value"] == 12345.0

    #setup
    model_settings["w1"] = 2.25
    #action
    common_params = create_common_parameters(
        Float64,
        precip_sources = model_settings["precip_sources"],
        precip_sinks = model_settings["precip_sinks"],
        Nd = model_settings["Nd"],
    )
    kid_params = create_kid_parameters(
        Float64,
        w1 = model_settings["w1"],
        t1 = model_settings["t1"],
        p0 = model_settings["p0"],
        qtot_flux_correction = model_settings["qtot_flux_correction"],
        r_dry = model_settings["r_dry"],
        std_dry = model_settings["std_dry"],
        κ = model_settings["κ"],
    )

    #test
    TT.@test common_params isa KD.Common.Parameters.CommonParameters
    TT.@test kid_params isa KD.K1DModel.Parameters.KinematicDriverParameters
    TT.@test kid_params.w1 == 2.25
end

TT.@testset "Run KiD" begin
    #setup
    u = [1e-4]
    u_names = ["cloud_liquid_water_specific_humidity_autoconversion_threshold"]
    model_settings = get_model_config()
    n_heights = [model_settings["n_elem"] for i in 1:2]
    n_times = length(model_settings["t_calib"])
    n_elem = model_settings["n_elem"]

    #action
    ode_sol, aux, precip = run_KiD(u, u_names, model_settings)

    #test
    TT.@test length(ode_sol) == n_times
    TT.@test length(parent(aux.thermo_variables.ρ_dry)) == model_settings["n_elem"]

    #action
    G = ODEsolution2Gvector(ode_sol, aux, precip, ["rlr", "rv"])

    #test
    TT.@test length(G) == sum(n_heights) * n_times

    #setup
    model_settings["filter"] = make_filter_props(
        [n_elem, n_elem],
        model_settings["t_calib"];
        apply = true,
        nz_per_filtered_cell = [2, 5],
        nt_per_filtered_cell = 5,
    )
    n_heights = model_settings["filter"]["nz_filtered"]

    #action
    ode_sol, aux = run_KiD(u, u_names, model_settings)
    G = ODEsolution2Gvector(ode_sol, aux, precip, ["rlr", "rv"], model_settings["filter"])

    #test
    TT.@test length(ode_sol) == (n_times - 1) * 5
    TT.@test length(G) == sum(n_heights) * (n_times - 1)

end

TT.@testset "Run KiD (multiple cases) and run dynamical model" begin
    #setup
    u = [1e-4]
    u_names = ["cloud_liquid_water_specific_humidity_autoconversion_threshold"]
    config = Dict("model" => get_model_config(), "observations" => get_observations_config())
    n_heights = config["model"]["n_elem"]
    n_times = length(config["model"]["t_calib"])
    n_vars = length(config["observations"]["data_names"])
    n_cases = length(config["observations"]["cases"])
    n_tot = n_heights * n_times * n_vars * n_cases
    case_numbers_1 = [1]
    case_numbers_tot = collect(1:n_cases)

    #action
    G_KiD_1 = run_KiD_multiple_cases(u, u_names, config, case_numbers_1)
    G_KiD = run_KiD_multiple_cases(u, u_names, config, case_numbers_tot)
    G_dyn_1 = run_dyn_model(u, u_names, config, case_numbers = case_numbers_1)
    G_dyn = run_dyn_model(u, u_names, config)

    #test
    TT.@test length(G_KiD_1) == n_heights * n_times * n_vars
    TT.@test length(G_KiD) == n_tot
    TT.@test length(G_dyn_1) == n_heights * n_times * n_vars
    TT.@test length(G_dyn) == n_tot
end
