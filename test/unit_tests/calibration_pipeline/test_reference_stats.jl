@testset "Test PCA computations" begin
    #setup
    y_mean = rand(10)
    covmat = convert(Array, Diagonal(shuffle(collect(range(0.1, 1.0, 10)))))
    allowed_var_loss = 0.1

    #action
    λ_pca, P_pca = KID.pca(covmat, allowed_var_loss)

    #test
    @test λ_pca == collect(range(0.3, 1.0, 8))
    @test maximum(P_pca, dims = 1) == maximum(P_pca, dims = 1) == ones(1, 8)

    #action
    y_pca, covmat_pca, P_pca = KID.obs_PCA(y_mean, covmat, allowed_var_loss)

    #action
    @test length(y_pca) == 8
    @test y_pca == P_pca' * y_mean
    @test covmat_pca == Diagonal(λ_pca)
    @test maximum(P_pca, dims = 1) == maximum(P_pca, dims = 1) == ones(1, 8)
end

@testset "Constructing reference statistics" begin
    #setup
    stats_config =
        Dict("perform_pca" => true, "variance_loss" => 0.1, "tikhonov_mode" => "absolute", "tikhonov_noise" => 0.0)
    y = rand(10, 50)
    covmat = convert(Array, Diagonal(shuffle(collect(range(0.1, 1.0, 10)))))
    truth = Observations.Observation(y, covmat, ["_"])

    #action
    ref_stats = KID.ReferenceStatistics(truth, stats_config)

    #test
    @test ref_stats.y_full == mean(y, dims = 2)[:]
    @test ref_stats.Γ_full == covmat
    @test length(ref_stats.y) == 8
    @test sum(ref_stats.Γ, dims = 1)[:] == diag(ref_stats.Γ)
    @test size(ref_stats.Γ) == (8, 8)
    @test size(ref_stats.P_pca) == (10, 8)

    #setup
    stats_config =
        Dict("perform_pca" => false, "variance_loss" => 0.1, "tikhonov_mode" => "absolute", "tikhonov_noise" => 0.01)

    #action
    ref_stats = KID.ReferenceStatistics(truth, stats_config)

    #test
    @test ref_stats.y_full == mean(y, dims = 2)[:]
    @test ref_stats.Γ_full == covmat
    @test ref_stats.y == ref_stats.y_full
    @test sum(ref_stats.Γ, dims = 1)[:] == diag(ref_stats.Γ)
    @test ref_stats.Γ == ref_stats.Γ_full .+ 0.01 * I(10)
    @test ref_stats.P_pca == I(10)
end
