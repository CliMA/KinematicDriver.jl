@testset "Get validation parameters" begin
    #setup
    config = get_config()
    vals = [v.mean for v in values(config["prior"]["parameters"])]
    α = 1 - config["observations"]["true_values_offset"]

    #action
    ϕ, θ = KID.params_validation(config)

    #test
    @test ϕ == α .* vals
    @test θ == 0.0

    #setup
    priors = KID.construct_priors(config["prior"]["parameters"])

    #action
    ϕ, θ = KID.params_validation(config, priors)

    #test
    @test ϕ == α .* vals
    @test length(θ) == length(vals)
end

@testset "Get validation samples" begin
    #setup
    config = get_config()
    n_heights = config["model"]["n_elem"]
    n_times = length(config["model"]["t_calib"])
    n_vars = length(config["observations"]["data_names"])
    n_cases = length(config["observations"]["cases"])
    n_tot = n_heights * n_times * n_vars * n_cases
    n_samples = config["observations"]["number_of_samples"]

    #action
    y = KID.get_validation_samples(config)

    #test
    @test size(y) == (n_tot, n_samples)
    @test norm(maximum(y, dims = 2) - mean(y, dims = 2)) < eps(Float64)

    #setup
    config["observations"]["scov_G_ratio"] = 0.2

    #action
    y = KID.get_validation_samples(config)

    #test
    @test size(y) == (n_tot, n_samples)
    @test sum((maximum(y, dims = 2) - mean(y, dims = 2)) .^ 2) > 1e-3
end

@testset "filter field" begin
    #setup
    a = [
        1 2 3 4 5 6
        7 8 9 1 2 3.0
        0 1 0 3 2 9
        -6 5 -3 4 2 1
    ]
    t_data = [1, 2, 3, 4, 5, 6, 7.0]
    t_calib = [1, 4, 6.0]
    z_data = [-4, -3, -2, -1, 0.0]
    z_calib = [-3, -2, 0.0]

    #action
    a_ = KID.filter_field(a, t_data, t_calib, z_data, z_calib)

    #test
    @test a_ == [8.0 1.5; -0.5 2.75]
end

@testset "Get observational data" begin
    #setup
    cases = get_observations_config()["cases"]
    dirs = [case.dir for case in cases]
    n_samples = 50
    n_heights = 50
    n_times = 31
    n_cases = length(dirs)
    generate_fake_pysdm_data(dirs; n_files = n_samples, z_max = 3000.0, n_z = 60, t_max = 240.0, n_t = 41)
    heights = collect(range(0.0, 3000, n_heights))
    times = collect(range(0.0, 240.0, n_times))
    variables = ["qlr", "rr"]
    n_variables = length(variables)
    n_single_case = n_heights * n_times * n_variables
    n_multiple_cases = n_single_case * n_cases

    #action
    result = KID.get_single_obs_field(dirs[1] * "output_1.nc", variables, heights, times)

    #test
    @test result isa Dict
    @test Set(keys(result)) == Set(variables)
    @test size(result["qlr"]) == size(result["rr"]) == (50, 31)
    @test_throws Exception KID.get_single_obs_field(
        dirs[1] * "output_1.nc",
        variables,
        heights,
        times;
        apply_filter = true,
    )


    #action
    result_vec = KID.get_single_obs_vector(dirs[1] * "output_1.nc", variables, heights, times)

    #test
    @test result_vec isa Vector
    @test length(result_vec) == n_single_case

    #action
    result_mat = KID.get_obs_matrix(dirs[1], variables, heights, times)

    #test
    @test result_mat isa Matrix
    @test size(result_mat) == (n_single_case, n_samples)

    #action
    result_tot = KID.get_obs_matrix(cases, variables, heights, times)

    #test
    @test result_tot isa Matrix
    @test size(result_tot) == (n_multiple_cases, n_samples)

    #setup
    n_heights = 30
    n_times = 21
    heights = collect(range(50.0, 2950, n_heights))
    times = collect(range(0.0, 240.0, n_times))
    n_single_case = n_heights * (n_times - 1) * n_variables
    n_multiple_cases = n_single_case * n_cases

    #action
    result_tot = KID.get_obs_matrix(cases, variables, heights, times; apply_filter = true)

    #test
    @test result_tot isa Matrix
    @test size(result_tot) == (n_multiple_cases, n_samples)

    rm_fake_pysdm_data(dirs)
end

@testset "normalize observations" begin
    # setup
    n_cases = 2
    n_heights = 2
    n_times = 2
    n_vars = 2
    n_samples = 20
    y = rand(n_cases * n_vars * n_heights * n_times, n_samples) .* 10.0
    norm_vecs = ones(n_cases, n_vars)

    # action
    ym = mean(y, dims = 2)
    _norm_vecs = [
        maximum([ym[1:2]; ym[5:6]]) maximum([ym[3:4]; ym[7:8]])
        maximum([ym[9:10]; ym[13:14]]) maximum([ym[11:12]; ym[15:16]])
    ]
    _norm_vecs_extended = reshape(repeat(_norm_vecs, inner = [2 2])', :, 1)
    _y = y ./ _norm_vecs_extended
    KID.normalize_y!(y, norm_vecs, n_vars, n_heights, n_times)

    #test
    @test norm_vecs == _norm_vecs
    @test y == _y
end

@testset "Get Observations" begin
    #setup
    config = get_config()
    n_heights = config["model"]["n_elem"]
    n_times = length(config["model"]["t_calib"])
    n_vars = length(config["observations"]["data_names"])
    n_cases = length(config["observations"]["cases"])
    n_tot = n_heights * n_times * n_vars * n_cases
    n_samples = config["observations"]["number_of_samples"]

    #action
    truth = KID.get_obs!(config)

    #test
    @test truth isa Observation
    @test length(truth.mean) == n_tot
    @test size(truth.obs_noise_cov) == (n_tot, n_tot)

    #setup
    config["observations"]["data_source"] = "file"
    dirs = [case.dir for case in config["observations"]["cases"]]
    generate_fake_pysdm_data(dirs; n_files = 50, z_max = 3000.0, n_z = 60, t_max = 240.0, n_t = 40)

    #action
    truth = KID.get_obs!(config)

    #test
    @test truth isa Observation
    @test length(truth.mean) == n_tot
    @test size(truth.obs_noise_cov) == (n_tot, n_tot)

    rm_fake_pysdm_data(dirs)
end
