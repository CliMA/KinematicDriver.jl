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
