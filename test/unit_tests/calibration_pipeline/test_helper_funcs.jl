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
    ext = KID.get_limits(u::Vector{Matrix{Float64}})
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
    filter = KID.make_filter_props(n_elem, t_calib; apply = true, nz_per_filtered_cell = 2, nt_per_filtered_cell = 2)

    #test
    @test filter["apply"] == true
    @test filter["nz_filtered"] == 5
    @test filter["nt_filtered"] == 2
    @test filter["saveat_t"] == [25.0, 75.0, 200.0, 400.0]
end

@testset "get numbers from config" begin
    #setup
    config = get_config()

    #action
    nums = KID.get_numbers_from_config(config)

    #test
    @test length(nums) == 4
    @test eltype(nums) == Int
end

@testset "make block diagonal matrix" begin
    #setup
    a = [1.0 2; 3 4]
    b = [5.0 6 7]

    #action
    c = KID.make_block_diagonal_matrix(a, b)

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
    obs = KID.get_obs!(config)
    ref_stats_list = KID.make_ref_stats_list(obs, config["statistics"], KID.get_numbers_from_config(config)...)
    ref_stats = KID.combine_ref_stats(ref_stats_list)

    #action
    me_1 = KID.compute_error_metrics(ϕ_values, ϕ_names, config, obs)
    me_2 = KID.compute_error_metrics(ϕ_values, ϕ_names, config, ref_stats)

    #test
    @test me_1.loss ≈ me_2.loss rtol = eps(Float64) * 10.0
    @test me_1.mse_m ≈ me_2.mse_m rtol = eps(Float64) * 10.0
    @test me_1.mse_s ≈ me_2.mse_s rtol = eps(Float64) * 10.0
end
