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
    rf_6 = "VarTimeScaleAcnv"
    st_1 = "CliMA_1M"
    st_2 = "Chen2022"

    #action
    moisture_eq, precip_n = KID.equation_types(eqm, pn, "_", "_")
    moisture_neq, precip_n = KID.equation_types(neqm, pn, "_", "_")
    moisture_eq, precip_0m = KID.equation_types(eqm, p0m, "_", "_")
    moisture_eq, precip_1m_1 = KID.equation_types(eqm, p1m, rf_1, st_1)
    moisture_eq, precip_1m_2 = KID.equation_types(eqm, p1m, rf_2, st_1)
    moisture_eq, precip_1m_3 = KID.equation_types(eqm, p1m, rf_3, st_1)
    moisture_eq, precip_1m_4 = KID.equation_types(eqm, p1m, rf_4, st_1)
    moisture_eq, precip_1m_5 = KID.equation_types(eqm, p1m, rf_5, st_1)
    moisture_eq, precip_1m_6 = KID.equation_types(eqm, p1m, rf_6, st_1)
    moisture_eq, precip_1m_7 = KID.equation_types(eqm, p1m, rf_6, st_2)
    #setup
    @test_throws Exception KID.equation_types("_", pn, rf_1, st_1)
    @test_throws Exception KID.equation_types(eqm, "_", rf_1, st_1)
    @test_throws Exception KID.equation_types(eqm, p1m, "_", st_1)
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

@testset "Terminal velocity" begin
    #setup
    u = Vector{Float64}(undef, 0)
    u_names = Vector{String}(undef, 0)

    n_h = 2
    n_t = 2
    n_v = 3
    n_c = 1

    model_settings = get_model_config(u_names)
    model_settings["n_elem"] = n_h
    model_settings["t_calib"] = collect(range(0, 60, n_t))

    n = n_h * n_t * n_v * n_c
    obs = rand(n)

    fixed_microphys_param_pairs = model_settings["fixed_microphys_param_pairs"]
    thermo_params = model_settings["thermo_params"]
    modal_nucleation_params = model_settings["modal_nucleation_params"]
    microphys_params = CM.Parameters.CloudMicrophysicsParameters(;
        fixed_microphys_param_pairs...,
        thermo_params,
        modal_nucleation_params,
    )
    ρ1 = obs[1:2]
    q_rai1 = obs[3:4]
    ρ2 = obs[7:8]
    q_rai2 = obs[9:10]

    #action
    vel = KID.compute_terminal_velocity(u, u_names, model_settings, obs)

    #test
    @test vel[5:6] ==
          CM.Microphysics1M.terminal_velocity.(
        microphys_params,
        CM.CommonTypes.RainType(),
        CM.CommonTypes.Blk1MVelType(),
        ρ1,
        q_rai1,
    )
    @test vel[11:12] ==
          CM.Microphysics1M.terminal_velocity.(
        microphys_params,
        CM.CommonTypes.RainType(),
        CM.CommonTypes.Blk1MVelType(),
        ρ2,
        q_rai2,
    )
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
    G = KID.ODEsolution2Gvector(ode_sol, aux, ["rlr", "rv"])

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
    G = KID.ODEsolution2Gvector(ode_sol, aux, ["rlr", "rv"], model_settings["filter"])

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

    #setup
    config["model"]["model"] = "terminal_velocity"
    config["observations"]["data_names"] = ["rho", "qr", "rain averaged terminal velocity"]
    n_tot = n_heights * n_times * 3 * n_cases
    config["statistics"] = get_stats_config()
    ref_stats = KID.combine_ref_stats(
        KID.make_ref_stats_list(rand(n_tot, 50), config["statistics"], n_cases, 3, n_heights, n_times),
    )

    #action
    G_vel = KID.run_dyn_model(u, u_names, config, RS = ref_stats)

    #test
    @test_throws Exception KID.run_dyn_model(u, u_names, config)
    @test length(G_vel) == n_tot
end
