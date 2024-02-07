@testset "creating KiD parameters object" begin
    #setup
    u = [1e-4, 12345.0]
    u_names = ["q_liq_threshold", "τ_acnv_rai"]
    model_settings = get_model_config()
    #action
    KID.update_parameters!(model_settings, u, u_names)
    #test
    @test model_settings["toml_dict"]["q_liq_threshold"]["value"] == 1e-4
    @test model_settings["toml_dict"]["τ_acnv_rai"]["value"] == 12345.0

    #setup
    model_settings["w1"] = 2.25
    #action
    kid_params = KID.create_kid_parameters(Float64, model_settings)
    #test
    @test kid_params isa KP.KinematicParameters
    @test kid_params.w1 == 2.25
end

@testset "Equation types" begin
    #setup
    toml_dict = get_model_config()["toml_dict"]
    eqm = "EquilibriumMoisture"
    neqm = "NonEquilibriumMoisture"
    pn = "NoPrecipitation"
    p0m = "Precipitation0M"
    p1m = "Precipitation1M"
    rf_1 = "CliMA_1M"
    rf_2 = "KK2000"
    rf_3 = "B1994"
    rf_4 = "LD2004"
    rf_5 = "TC1980"
    rf_6 = "VarTimeScaleAcnv"
    st_1 = "CliMA_1M"
    st_2 = "Chen2022"

    #action
    moisture_eq, precip_n = KID.equation_types(eqm, pn, "_", "_", toml_dict)
    moisture_neq, precip_n = KID.equation_types(neqm, pn, "_", "_", toml_dict)
    moisture_eq, precip_0m = KID.equation_types(eqm, p0m, "_", "_", toml_dict)
    moisture_eq, precip_1m_1 = KID.equation_types(eqm, p1m, rf_1, st_1, toml_dict)
    moisture_eq, precip_1m_2 = KID.equation_types(eqm, p1m, rf_2, st_1, toml_dict)
    moisture_eq, precip_1m_3 = KID.equation_types(eqm, p1m, rf_3, st_1, toml_dict)
    moisture_eq, precip_1m_4 = KID.equation_types(eqm, p1m, rf_4, st_1, toml_dict)
    moisture_eq, precip_1m_5 = KID.equation_types(eqm, p1m, rf_5, st_1, toml_dict)
    moisture_eq, precip_1m_6 = KID.equation_types(eqm, p1m, rf_6, st_1, toml_dict)
    moisture_eq, precip_1m_7 = KID.equation_types(eqm, p1m, rf_6, st_2, toml_dict)
    #setup
    @test_throws Exception KID.equation_types("_", pn, rf_1, st_1, toml_dict)
    @test_throws Exception KID.equation_types(eqm, "_", rf_1, st_1, toml_dict)
    @test_throws Exception KID.equation_types(eqm, p1m, "_", st_1, toml_dict)
    @test moisture_eq isa KID.EquilibriumMoisture_ρdTq
    @test moisture_neq isa KID.NonEquilibriumMoisture_ρdTq
    @test precip_n isa KID.NoPrecipitation
    @test precip_0m isa KID.Precipitation0M
    @test precip_1m_1 isa KID.Precipitation1M
    @test precip_1m_2 isa KID.Precipitation1M
    @test precip_1m_3 isa KID.Precipitation1M
    @test precip_1m_4 isa KID.Precipitation1M
    @test precip_1m_5 isa KID.Precipitation1M
    @test precip_1m_6 isa KID.Precipitation1M
    @test precip_1m_7 isa KID.Precipitation1M
end

@testset "Run KiD" begin
    #setup
    u = [1e-4]
    u_names = ["q_liq_threshold"]
    model_settings = get_model_config()
    n_heights = model_settings["n_elem"]
    n_times = length(model_settings["t_calib"])

    #action
    ode_sol, aux, precip = KID.run_KiD(u, u_names, model_settings)

    #test
    @test length(ode_sol) == n_times
    @test length(parent(aux.moisture_variables.ρ_dry)) == n_heights

    #action
    G = KID.ODEsolution2Gvector(ode_sol, aux, precip, ["rlr", "rv"])

    #test
    @test length(G) == n_heights * n_times * 2

    #setup
    model_settings["filter"] = KID.make_filter_props(
        model_settings["n_elem"],
        model_settings["t_calib"];
        apply = true,
        nz_per_filtered_cell = 2,
        nt_per_filtered_cell = 5,
    )

    #action
    ode_sol, aux = KID.run_KiD(u, u_names, model_settings)
    G = KID.ODEsolution2Gvector(ode_sol, aux, precip, ["rlr", "rv"], model_settings["filter"])

    #test
    @test_throws Exception KID.get_variable_data_from_ODE(ode_sol[1], aux, precip, "qt_")
    for var in ["ql", "rl", "rho", "rain averaged terminal velocity", "rainrate"]
        @test length(KID.get_variable_data_from_ODE(ode_sol[1], aux, precip, "rain averaged terminal velocity")) ==
              n_heights
    end
    @test length(ode_sol) == (n_times - 1) * 5
    @test length(G) == n_heights / 2 * (n_times - 1) * 2

end

@testset "Run KiD (multiple cases) and run dynamical model" begin
    #setup
    u = [1e-4]
    u_names = ["q_liq_threshold"]
    config = Dict("model" => get_model_config(), "observations" => get_observations_config())
    n_heights = config["model"]["n_elem"]
    n_times = length(config["model"]["t_calib"])
    n_vars = length(config["observations"]["data_names"])
    n_cases = length(config["observations"]["cases"])
    n_tot = n_heights * n_times * n_vars * n_cases
    case_numbers_1 = [1]
    case_numbers_tot = collect(1:n_cases)

    #action
    G_KiD_1 = KID.run_KiD_multiple_cases(u, u_names, config, case_numbers_1)
    G_KiD = KID.run_KiD_multiple_cases(u, u_names, config, case_numbers_tot)
    G_dyn_1 = KID.run_dyn_model(u, u_names, config, case_numbers = case_numbers_1)
    G_dyn = KID.run_dyn_model(u, u_names, config)

    #test
    @test length(G_KiD_1) == n_heights * n_times * n_vars
    @test length(G_KiD) == n_tot
    @test length(G_dyn_1) == n_heights * n_times * n_vars
    @test length(G_dyn) == n_tot
end
