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
