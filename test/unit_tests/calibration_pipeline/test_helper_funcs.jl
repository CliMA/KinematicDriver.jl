@testset "Find limits" begin
    #setup
    n_vars = 10
    n_ens = 20
    n_iter = 30
    ensemble_matrix = rand(n_vars, n_ens)
    u = [ensemble_matrix .+ (i - 1) for i in 1:n_iter]
    minima = minimum(ensemble_matrix, dims = 2)
    maxima = maximum(ensemble_matrix, dims = 2) .+ (n_iter - 1)

    #action
    ext = KCP.get_limits(u::Vector{Matrix{Float64}})
    ext_min = [e[1] for e in ext]
    ext_max = [e[2] for e in ext]

    #test
    @test ext isa Vector{Tuple{Float64, Float64}}
    @test length(ext) == n_vars
    @test norm(ext_min .- minima) < eps(Float64)
    @test norm(ext_max .- maxima) < eps(Float64)
end

@testset "make filter" begin
    #setup
    n_elem = 10
    t_calib = [0, 100, 500]

    #action
    filter =
        KCP.make_filter_props([n_elem], t_calib; apply = true, nz_per_filtered_cell = [2], nt_per_filtered_cell = 2)

    #test
    @test filter["apply"] == true
    @test filter["nz_filtered"] == [5]
    @test filter["nt_filtered"] == 2
    @test filter["saveat_t"] == [25.0, 75.0, 200.0, 400.0]
end

@testset "get numbers from config" begin
    #setup
    config = get_config()

    #action
    nums = KCP.get_numbers_from_config(config)

    #test
    @test length(nums) == 3
    @test typeof(nums) == @NamedTuple{n_cases::Int, n_heights::Vector{Int}, n_times::Int}
    @test length(nums.n_heights) == length(config["observations"]["data_names"])
end

@testset "get single case data from vector" begin
    #setup
    config = get_config()
    config["observations"]["data_names"] = ["rho", "ql", "qr", "rainrate"]
    config["model"]["filter"] =
        KCP.make_filter_props(config["model"]["n_elem"] .* ones(Int, 4), config["model"]["t_calib"])
    (n_c, n_z, n_t) = KCP.get_numbers_from_config(config)
    n_single_case = sum(n_z) * n_t
    vec = rand(n_c * n_single_case)

    #action
    single_case_vec = KCP.get_case_i_vec(vec, 1, n_single_case)
    fields = KCP.get_single_case_fields(single_case_vec, n_z, n_t)

    #test
    length(single_case_vec) == n_single_case
    for i in 1:length(n_z)
        size(fields[i]) == (n_z[i], n_t)
    end
end

@testset "make block diagonal matrix" begin
    #setup
    a = [1.0 2; 3 4]
    b = [5.0 6 7]

    #action
    c = KCP.make_block_diagonal_matrix(a, b)

    #test
    @test size(c) == (3, 5)
    @test c[1:2, 1:2] == a
    @test c[3:3, 3:5] == b
    @test c[1:2, 3:5] == zeros((2, 3))
    @test c[3:3, 1:2] == zeros((1, 2))
end

@testset "compute error metrics" begin
    #setup
    config = get_config()
    ϕ_names = collect(keys(config["prior"]["parameters"]))
    ϕ_values = collect([v.mean for v in values(config["prior"]["parameters"])])
    obs = KCP.get_obs!(config)
    ref_stats_list = KCP.make_ref_stats_list(obs, config["statistics"], KCP.get_numbers_from_config(config)...)
    ref_stats = KCP.combine_ref_stats(ref_stats_list)

    #action
    me_1 = KCP.compute_error_metrics(ϕ_values, ϕ_names, config, obs)
    me_2 = KCP.compute_error_metrics(ϕ_values, ϕ_names, config, ref_stats)

    #test
    @test me_1.loss ≈ me_2.loss rtol = eps(Float64) * 10.0
    @test me_1.mse_m ≈ me_2.mse_m rtol = eps(Float64) * 10.0
    @test me_1.mse_s ≈ me_2.mse_s rtol = eps(Float64) * 10.0
end
