@testset "Calibrate by Ensemble Kalman Processes - Get results" begin
    #setup
    config = get_config()
    priors = KID.construct_priors(config["prior"]["parameters"])
    truth = KID.get_obs!(config)
    ref_stats = KID.ReferenceStatistics(truth, config["statistics"])
    min_iter = config["process"]["n_iter_min"]
    max_iter = config["process"]["n_iter_max"]
    n_ensemble = config["process"]["n_ens"]
    n_vars = length(keys(config["prior"]["parameters"]))
    α = (1 - config["observations"]["true_values_offset"])
    u_true = [v.mean for v in values(config["prior"]["parameters"])] .* α

    #action
    res = KID.calibrate(KID.EKPStyle(), priors, config, ref_stats, verbose = false)
    θ = KID.get_u(res)
    u = KID.get_results(res, priors)

    #test
    @test res isa EnsembleKalmanProcess
    @test min_iter <= length(θ) - 1 <= max_iter
    @test size(θ[1]) == (n_vars, n_ensemble)
    @test u isa Vector
    @test norm((u .- u_true) ./ u_true) < 3e-1
end

@testset "Calibrate by Optim - Get results" begin
    #setup
    config = get_config()
    priors = KID.construct_priors(config["prior"]["parameters"])
    truth = KID.get_obs!(config)
    ref_stats = KID.ReferenceStatistics(truth, config["statistics"])
    u_names = collect(keys(config["prior"]["parameters"]))
    α = (1 - config["observations"]["true_values_offset"])
    u_true = [v.mean for v in values(config["prior"]["parameters"])] .* α

    #action
    res = KID.calibrate(KID.OptimStyle(), priors, config, ref_stats, verbose = false)
    f = KID.make_optim_loss_function(u_names, priors, config, ref_stats)
    u = KID.get_results(res, priors)

    #test
    @test res isa Optim.MultivariateOptimizationResults
    @test f isa Function
    @test u isa Vector
    @test norm((u .- u_true) ./ u_true) < 1e-2
end

@testset "Calibrate function for abstract optimization style" begin
    #setup
    config = get_config()
    priors = KID.construct_priors(config["prior"]["parameters"])
    truth = KID.get_obs!(config)
    ref_stats = KID.ReferenceStatistics(truth, config["statistics"])

    #test
    @test_throws Exception KID.calibrate(DummyStyle(), priors, config, ref_stats, verbose = false)
end

@testset "Compute loss" begin
    #setup
    config = get_config()
    truth = KID.get_obs!(config)
    ref_stats = KID.ReferenceStatistics(truth, config["statistics"])
    u_names = collect(keys(config["prior"]["parameters"]))
    u = [v.mean for v in values(config["prior"]["parameters"])]
    α = (1 - config["observations"]["true_values_offset"])
    u_true = u .* α

    #action
    loss_u = KID.compute_loss(u, u_names, config, ref_stats)
    loss_u_true = KID.compute_loss(u_true, u_names, config, ref_stats)

    #test
    @test 0.0 < loss_u
    @test 0.0 <= loss_u_true < eps(Float64)
end
