TT.@testset "Calibrate by Ensemble Kalman Processes and Get results" begin
    #setup
    config = get_config()
    priors = construct_priors(config["prior"]["parameters"])
    obs = get_obs!(config)
    ref_stats_list = make_ref_stats_list(obs, config["statistics"], get_numbers_from_config(config)...)
    n_iter = config["process"]["n_iter"]
    n_cases = length(ref_stats_list)
    batch_size = config["process"]["batch_size"]
    n_batches = floor(Int, n_cases / batch_size)
    n_ensemble = config["process"]["n_ens"]
    n_vars = length(keys(config["prior"]["parameters"]))
    α = (1 - config["observations"]["true_values_offset"])
    u_true = [v.mean for v in values(config["prior"]["parameters"])] .* α

    #action
    res = calibrate(EKPStyle(), priors, config, ref_stats_list, verbose = 0)
    u = get_results(res, priors)

    #test
    TT.@test res isa Vector
    TT.@test length(res) == n_iter * n_batches + 1
    TT.@test size(res[1]) == (n_vars, n_ensemble)
    TT.@test u.ϕ_optim isa Vector
    TT.@test u.cov_optim isa Matrix
    TT.@test LA.norm((u.ϕ_optim .- u_true) ./ u_true) < 0.1

    #setup
    config["process"]["EKP_method"] = "UKI"
    config["process"]["Δt"] = n_cases / batch_size

    #action
    res = calibrate(EKPStyle(), priors, config, ref_stats_list, verbose = 0)
    u = get_results(res, priors)

    #test
    TT.@test res isa Vector
    TT.@test length(res) == n_iter * n_batches + 1
    TT.@test size(res[1]) == (n_vars, 2 * n_vars + 1)
    TT.@test u.ϕ_optim isa Vector
    TT.@test u.cov_optim isa Matrix
    TT.@test LA.norm((u.ϕ_optim .- u_true) ./ u_true) < 0.13 #TODO

    #setup
    config["process"]["EKP_method"] = "UKP_"

    #test
    TT.@test_throws Exception calibrate(EKPStyle(), priors, config, ref_stats_list)
end

TT.@testset "EKP with mini-batching" begin
    #setup
    config = get_config()
    config["prior"]["parameters"] = Dict(
        "a" => (mean = 1.0, var = 0.5, lbound = -Inf, ubound = Inf),
        "b" => (mean = 5.0, var = 0.5, lbound = -Inf, ubound = Inf),
    )
    config["observations"]["data_names"] = ["test data"]
    config["model"]["filter"] =
        make_filter_props(config["model"]["n_elem"] .* ones(Int, 1), config["model"]["t_calib"])
    config["observations"]["cases"] = [
        (power = 1.0,),
        (power = 2.0,),
        (power = 3.0,),
        (power = 4.0,),
        (power = 5.0,),
        (power = 6.0,),
        (power = 7.0,),
        (power = 8.0,),
        (power = 9.0,),
        (power = 10.0,),
    ]
    config["model"]["model"] = "test_model"
    α = (1 - config["observations"]["true_values_offset"])
    u_true = [v.mean for v in values(config["prior"]["parameters"])] .* α

    priors = construct_priors(config["prior"]["parameters"])
    obs = get_obs!(config)
    ref_stats_list = make_ref_stats_list(obs, config["statistics"], get_numbers_from_config(config)...)

    #action
    res = calibrate(EKPStyle(), priors, config, ref_stats_list, verbose = 0)
    u = get_results(res, priors)

    #test
    TT.@test LA.norm((u.ϕ_optim .- u_true) ./ u_true) < 0.1

    #setup
    config["process"]["EKP_method"] = "UKI"
    config["process"]["Δt"] = 10.0

    #action
    res = calibrate(EKPStyle(), priors, config, ref_stats_list, verbose = 0)
    u = get_results(res, priors)

    #test
    TT.@test LA.norm((u.ϕ_optim .- u_true) ./ u_true) < 0.1
end

TT.@testset "Calibrate by Optim and Get results" begin
    #setup
    config = get_config()
    priors = construct_priors(config["prior"]["parameters"])
    obs = get_obs!(config)
    ref_stats_list = make_ref_stats_list(obs, config["statistics"], get_numbers_from_config(config)...)
    ref_stats = combine_ref_stats(ref_stats_list)
    u_names = collect(keys(config["prior"]["parameters"]))
    α = (1 - config["observations"]["true_values_offset"])
    u_true = [v.mean for v in values(config["prior"]["parameters"])] .* α

    #action
    res = calibrate(OptimStyle(), priors, config, ref_stats, verbose = 0)
    f = make_optim_loss_function(u_names, priors, config, ref_stats)
    u = get_results(res, priors)

    #test
    TT.@test res isa Optim.MultivariateOptimizationResults
    TT.@test f isa Function
    TT.@test u.ϕ_optim isa Vector
    TT.@test u.cov_optim isa Matrix
    TT.@test LA.norm((u.ϕ_optim .- u_true) ./ u_true) < 0.05
end

TT.@testset "Calibrate function for abstract optimization style" begin
    #setup
    config = get_config()
    priors = construct_priors(config["prior"]["parameters"])
    obs = get_obs!(config)
    ref_stats_list = make_ref_stats_list(obs, config["statistics"], get_numbers_from_config(config)...)
    ref_stats = combine_ref_stats(ref_stats_list)

    #test
    TT.@test_throws Exception calibrate(DummyStyle(), priors, config, ref_stats, verbose = 0)
end

TT.@testset "Compute loss" begin
    #setup
    config = get_config()
    config["observations"]["scov_G_ratio"] = 0.0
    obs = get_obs!(config)
    ref_stats_list = make_ref_stats_list(obs, config["statistics"], get_numbers_from_config(config)...)
    ref_stats = combine_ref_stats(ref_stats_list)
    u_names = collect(keys(config["prior"]["parameters"]))
    u = [v.mean for v in values(config["prior"]["parameters"])]
    α = (1 - config["observations"]["true_values_offset"])
    u_true = u .* α

    #action
    loss_u = compute_loss(u, u_names, config, ref_stats)
    loss_u_true = compute_loss(u_true, u_names, config, ref_stats)

    #test
    TT.@test eps(Float64) * 10.0 < loss_u
    TT.@test 0.0 <= loss_u_true < eps(Float64) * 10.0
end
