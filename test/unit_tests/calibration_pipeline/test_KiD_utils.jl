@testset "creating KiD parameters object" begin
    #setup
    u = [1e-4, 12345.0]
    u_names = ["q_liq_threshold", "τ_acnv_rai"]
    #action
    params_calib = KID.create_param_dict(u, u_names)
    #test
    @test length(keys(params_calib)) == 2
    @test params_calib["q_liq_threshold"] == 1e-4
    @test params_calib["τ_acnv_rai"] == 12345.0

    #setup
    model_settings = get_model_config(u_names)
    model_settings["w1"] = 2.25
    #action
    params = KID.create_parameter_set(Float64, model_settings, params_calib)
    #test
    @test params isa KP.KinematicParameters
    @test CM.Parameters.q_liq_threshold(KP.microphysics_params(params)) == 1e-4
    @test CM.Parameters.τ_acnv_rai(KP.microphysics_params(params)) == 12345.0
    @test KP.w1(params) == 2.25
end

@testset "Equation types" begin
    #setup
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

    #action
    moisture_eq, precip_n = KID.equation_types(eqm, pn, "_")
    moisture_neq, precip_n = KID.equation_types(neqm, pn, "_")
    moisture_eq, precip_0m = KID.equation_types(eqm, p0m, "_")
    moisture_eq, precip_1m_1 = KID.equation_types(eqm, p1m, rf_1)
    moisture_eq, precip_1m_2 = KID.equation_types(eqm, p1m, rf_2)
    moisture_eq, precip_1m_3 = KID.equation_types(eqm, p1m, rf_3)
    moisture_eq, precip_1m_4 = KID.equation_types(eqm, p1m, rf_4)
    moisture_eq, precip_1m_5 = KID.equation_types(eqm, p1m, rf_5)
    #setup
    @test_throws Exception KID.equation_types("_", pn, rf_1)
    @test_throws Exception KID.equation_types(eqm, "_", rf_1)
    @test_throws Exception KID.equation_types(eqm, p1m, "_")
    @test moisture_eq isa KID.EquilibriumMoisture_ρdTq
    @test moisture_neq isa KID.NonEquilibriumMoisture_ρdTq
    @test precip_n isa KID.NoPrecipitation
    @test precip_0m isa KID.Precipitation0M
    @test precip_1m_1 isa KID.Precipitation1M
    @test precip_1m_2 isa KID.Precipitation1M
    @test precip_1m_3 isa KID.Precipitation1M
    @test precip_1m_4 isa KID.Precipitation1M
    @test precip_1m_5 isa KID.Precipitation1M
end

@testset "Terminal velocity" begin
    #setup
    u = Vector{Float64}(undef, 0)
    u_names = Vector{String}(undef, 0)

    n_h = 2
    n_t = 2
    n_v = 3
    n_c = 1

    config = Dict("model" => get_model_config(u_names), "observations" => Dict("ynorm" => [1.0 1e-3 1.0]))
    config["model"]["n_elem"] = n_h
    config["model"]["t_calib"] = collect(range(0, 60, n_t))

    n = n_h * n_t * n_v * n_c
    obs = rand(n)

    fixed_microphys_param_pairs = config["model"]["fixed_microphys_param_pairs"]
    thermo_params = config["model"]["thermo_params"]
    microphys_params = CM.Parameters.CloudMicrophysicsParameters(; fixed_microphys_param_pairs..., thermo_params)
    ρ1 = obs[1:2] * config["observations"]["ynorm"][1, 1]
    q_rai1 = obs[3:4] * config["observations"]["ynorm"][1, 2]
    ρ2 = obs[7:8] * config["observations"]["ynorm"][1, 1]
    q_rai2 = obs[9:10] * config["observations"]["ynorm"][1, 2]

    #action
    vel = KID.compute_terminal_velocity(u, u_names, config, obs)

    #test
    @test vel[5:6] ==
          CM.Microphysics1M.terminal_velocity.(microphys_params, CM.CommonTypes.RainType(), ρ1, q_rai1) /
          config["observations"]["ynorm"][1, 3]
    @test vel[11:12] ==
          CM.Microphysics1M.terminal_velocity.(microphys_params, CM.CommonTypes.RainType(), ρ2, q_rai2) /
          config["observations"]["ynorm"][1, 3]
end

@testset "Get variable data from ODE" begin
    #setup
    ρ = [1.0, 1.0]
    u = (; ρq_tot = [0.01, 0.005], ρq_liq = [0.001, 0.0005], ρq_rai = [0.0001, 0.00005])

    #test
    @test_throws Exception KID.get_variable_data_from_ODE(u, ρ, "qt_")
    @test KID.get_variable_data_from_ODE(u, ρ, "ql") == [0.001, 0.0005]
    @test KID.get_variable_data_from_ODE(u, ρ, "rr") ≈ [0.00010101, 0.00005025] atol = 1e-8
    @test KID.get_variable_data_from_ODE(u, ρ, "rho") == ρ
    @test length(KID.get_variable_data_from_ODE(u, ρ, "rain averaged terminal velocity")) == length(ρ)
end

@testset "Run KiD" begin
    #setup
    u = [1e-4]
    u_names = ["q_liq_threshold"]
    model_settings = get_model_config(u_names)
    n_heights = model_settings["n_elem"]
    n_times = length(model_settings["t_calib"])

    #action
    ode_sol, aux = KID.run_KiD(u, u_names, model_settings)

    #test
    @test length(ode_sol) == n_times
    @test length(parent(aux.moisture_variables.ρ_dry)) == n_heights

    #action
    G = KID.ODEsolution2Gvector(ode_sol, aux, ["rlr", "rv"], [0.01, 1.0])

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
    G = KID.ODEsolution2Gvector(ode_sol, aux, ["rlr", "rv"], [0.01, 1.0], model_settings["filter"])

    #test
    @test length(ode_sol) == (n_times - 1) * 5
    @test length(G) == n_heights / 2 * (n_times - 1) * 2

end

@testset "Run KiD (multiple cases) and run dynamical model" begin
    #setup
    u = [1e-4]
    u_names = ["q_liq_threshold"]
    config = Dict("model" => get_model_config(u_names), "observations" => get_observations_config())
    n_heights = config["model"]["n_elem"]
    n_times = length(config["model"]["t_calib"])
    n_vars = length(config["observations"]["data_names"])
    n_cases = length(config["observations"]["cases"])
    n_tot = n_heights * n_times * n_vars * n_cases

    #action
    G_KiD = KID.run_KiD_multiple_cases(u, u_names, config)
    G_dyn = KID.run_dyn_model(u, u_names, config)

    #test
    @test length(G_KiD) == n_tot
    @test length(G_dyn) == n_tot

    #setup
    config["model"]["model"] = "terminal_velocity"
    config["observations"]["data_names"] = ["rho", "qr", "rain averaged terminal velocity"]
    config["observations"]["ynorm"] = ones(n_cases, 3)
    n_tot = n_heights * n_times * 3 * n_cases
    config["statistics"] = get_stats_config()
    ref_stats = KID.ReferenceStatistics(
        Observations.Observation(rand(n_tot, 50), rand(n_tot, n_tot), ["_"]),
        config["statistics"],
    )

    #action
    G_vel = KID.run_dyn_model(u, u_names, config, RS = ref_stats)

    #test
    @test_throws Exception KID.run_dyn_model(u, u_names, config)
    @test length(G_vel) == n_tot
end
