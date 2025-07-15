TT.@testset "Get validation parameters" begin
    #setup
    config = get_config()
    vals = [v.mean for v in values(config["prior"]["parameters"])]
    α = 1 - config["observations"]["true_values_offset"]

    #action
    ϕ, θ = params_validation(config)

    #test
    TT.@test ϕ == α .* vals
    TT.@test θ == 0.0

    #setup
    priors = construct_priors(config["prior"]["parameters"])

    #action
    ϕ, θ = params_validation(config, priors)

    #test
    TT.@test ϕ == α .* vals
    TT.@test length(θ) == length(vals)
end

TT.@testset "Get validation samples" begin
    #setup
    config = get_config()
    n_heights = config["model"]["n_elem"]
    n_times = length(config["model"]["t_calib"])
    n_vars = length(config["observations"]["data_names"])
    n_cases = length(config["observations"]["cases"])
    n_tot = n_heights * n_times * n_vars * n_cases
    n_samples = config["observations"]["number_of_samples"]
    config["observations"]["scov_G_ratio"] = 0.0

    #action
    y = get_validation_samples(config)

    #test
    TT.@test size(y) == (n_tot, n_samples)
    TT.@test LA.norm(maximum(y, dims = 2) - Statistics.mean(y, dims = 2)) < eps(Float64) * 10.0

    #setup
    config["observations"]["scov_G_ratio"] = 0.2

    #action
    y = get_validation_samples(config)

    #test
    TT.@test size(y) == (n_tot, n_samples)
    TT.@test sum((maximum(y, dims = 2) - Statistics.mean(y, dims = 2)) .^ 2) > eps(Float64) * 10.0
end

TT.@testset "filter field" begin
    #setup
    a = [
        1 2 3 4 5 6
        7 8 9 1 2 3.0
        0 1 0 3 2 9
        -6 5 -3 4 2 1
    ]
    t_data = [1, 2, 3, 4, 5, 6, 7.0]
    z_data = [-4, -3, -2, -1, 0.0]

    #test
    TT.@test_throws Exception filter_field(a, t_data, [1.5, 1.25], z_data, [-1.25, -0.5])
    TT.@test_throws Exception filter_field(a, t_data, [1.5, 1.5], z_data, [-1.25, -0.5])
    TT.@test_throws Exception filter_field(a, t_data, [1.0, 1.0], z_data, [-1.25, -0.5])
    TT.@test_throws Exception filter_field(a, t_data, [1.5, 3.25], z_data, [-3.0, -3.0])
    TT.@test filter_field(a, t_data, [1.5, 3.25], z_data, [-3.25, -0.5]) ≈ [3.3636363636363638] atol =
        eps(Float64) * 10.0
    TT.@test filter_field(a, t_data, [1, 4, 6.0], z_data, [-3, -2, 0.0]) ≈ [8.0 1.5; -0.5 2.75] atol =
        eps(Float64) * 10.0
end

TT.@testset "Get observational data" begin
    #setup
    cases = get_observations_config()["cases"]
    dirs = [case.dir for case in cases]
    n_samples = 50
    n_heights = [50, 50]
    n_times = 31
    n_cases = length(dirs)
    generate_fake_pysdm_data(dirs; n_files = n_samples, z_max = 3000.0, n_z = 60, t_max = 240.0, n_t = 41)
    heights = collect.(range.(0.0, 3000, n_heights))
    times = collect(range(0.0, 240.0, n_times))
    variables = ["qlr", "rr"]
    n_variables = length(variables)
    n_single_case = sum(n_heights) * n_times
    n_multiple_cases = n_single_case * n_cases

    #action
    result = get_single_obs_field(dirs[1] * "output_1.nc", variables, heights, times)

    #test
    TT.@test result isa Dict
    TT.@test Set(keys(result)) == Set(variables)
    TT.@test size(result["qlr"]) == size(result["rr"]) == (50, 31)
    TT.@test_throws Exception get_single_obs_field(
        dirs[1] * "output_1.nc",
        variables,
        heights,
        times;
        apply_filter = true,
    )


    #action
    result_vec = get_single_obs_vector(dirs[1] * "output_1.nc", variables, heights, times)

    #test
    TT.@test result_vec isa Vector
    TT.@test length(result_vec) == n_single_case

    #action
    result_mat = get_obs_matrix(dirs[1], variables, heights, times)

    #test
    TT.@test result_mat isa Matrix
    TT.@test size(result_mat) == (n_single_case, n_samples)

    #action
    result_tot = get_obs_matrix(cases, variables, heights, times)

    #test
    TT.@test result_tot isa Matrix
    TT.@test size(result_tot) == (n_multiple_cases, n_samples)

    #setup
    n_heights = [30, 30]
    n_times = 21
    heights = collect.(range.(50.0, 2950, n_heights))
    times = collect(range(0.0, 240.0, n_times))
    n_single_case = sum(n_heights) * (n_times - 1)
    n_multiple_cases = n_single_case * n_cases

    #action
    result_tot = get_obs_matrix(cases, variables, heights, times; apply_filter = true)

    #test
    TT.@test result_tot isa Matrix
    TT.@test size(result_tot) == (n_multiple_cases, n_samples)

    rm_fake_pysdm_data(dirs)
end

TT.@testset "Get Observations" begin
    #setup
    config = get_config()
    n_heights = [config["model"]["n_elem"] for i in 1:2]
    n_times = length(config["model"]["t_calib"])
    n_vars = length(config["observations"]["data_names"])
    n_cases = length(config["observations"]["cases"])
    n_tot = sum(n_heights) * n_times * n_cases
    n_samples = config["observations"]["number_of_samples"]

    #action
    obs = get_obs!(config)

    #test
    TT.@test obs isa Matrix
    TT.@test size(obs) == (n_tot, n_samples)

    #setup
    config["observations"]["data_source"] = "file"
    dirs = [case.dir for case in config["observations"]["cases"]]
    generate_fake_pysdm_data(dirs; n_files = 50, z_max = 3000.0, n_z = 60, t_max = 240.0, n_t = 40)

    #action
    truth = get_obs!(config)

    #test
    TT.@test obs isa Matrix
    TT.@test size(obs) == (n_tot, n_samples)

    rm_fake_pysdm_data(dirs)
end
