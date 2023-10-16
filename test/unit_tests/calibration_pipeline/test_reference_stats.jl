@testset "Test PCA computations" begin
    #setup
    y_mean = rand(10)
    covmat = convert(Array, Diagonal(shuffle(collect(range(0.1, 1.0, 10)))))
    allowed_var_loss = 0.1

    #action
    λ_pca, P_pca = KCP.pca(covmat, allowed_var_loss)

    #test
    @test λ_pca == collect(range(0.3, 1.0, 8))
    @test maximum(P_pca, dims = 1) == maximum(P_pca, dims = 1) == ones(1, 8)

    #action
    y_pca, covmat_pca, P_pca = KCP.obs_PCA(y_mean, covmat, allowed_var_loss)

    #action
    @test length(y_pca) == 8
    @test y_pca == P_pca' * y_mean
    @test covmat_pca == Diagonal(λ_pca)
    @test maximum(P_pca, dims = 1) == maximum(P_pca, dims = 1) == ones(1, 8)
end

@testset "Constructing reference statistics" begin
    #setup
    stats_config = Dict(
        "normalization" => "mean_normalized",
        "perform_pca" => true,
        "variance_loss" => 0.1,
        "tikhonov_mode" => "absolute",
        "tikhonov_noise" => 0.0,
    )
    n_cases = 1
    n_variables = 2
    n_heights = 3
    n_times = 2
    n_tot = n_cases * n_variables * n_heights * n_times
    n_samples = 50
    case_number = 1
    obs = rand(n_tot, n_samples)
    covmat = cov(obs, dims = 2)

    #action
    ref_stats = KCP.ReferenceStatistics(obs, stats_config, case_number, n_variables, n_heights, n_times)

    #test
    @test length(ref_stats.y_full) == n_tot
    @test size(ref_stats.Γ_full) == (n_tot, n_tot)
    @test length(ref_stats.y) < n_tot
    @test sum(ref_stats.Γ, dims = 1)[:] == diag(ref_stats.Γ)
    @test size(ref_stats.Γ)[1] == size(ref_stats.Γ)[2] < n_tot
    @test size(ref_stats.P_pca)[1] == n_tot
    @test size(ref_stats.P_pca)[2] < n_tot
    @test ref_stats.obs_mean ≈ mean(obs, dims = 2) rtol = eps(Float64) * 10.0
    @test size(ref_stats.y_norm) == (n_cases, n_variables)
    @test ref_stats.y_norm == ref_stats.var_mean_max
    @test size(ref_stats.var_std_max) == (n_cases, n_variables)
    @test ref_stats.centered == false
    @test ref_stats.case_numbers == [case_number]
    @test ref_stats.n_variables == n_variables
    @test ref_stats.n_heights == n_heights
    @test ref_stats.n_times == n_times

    #setup
    stats_config = Dict(
        "normalization" => "mean_normalized",
        "perform_pca" => false,
        "variance_loss" => 0.1,
        "tikhonov_mode" => "absolute",
        "tikhonov_noise" => 0.01,
    )

    #action
    ref_stats = KCP.ReferenceStatistics(obs, stats_config, case_number, n_variables, n_heights, n_times)

    #test
    @test length(ref_stats.y_full) == n_tot
    @test size(ref_stats.Γ_full) == (n_tot, n_tot)
    @test ref_stats.y == ref_stats.y_full
    @test ref_stats.Γ == ref_stats.Γ_full .+ 0.01 * I(n_tot)
    @test ref_stats.P_pca == I(n_tot)
    @test ref_stats.obs_mean ≈ mean(obs, dims = 2) rtol = eps(Float64) * 10.0
    @test size(ref_stats.y_norm) == (n_cases, n_variables)
    @test ref_stats.y_norm == ref_stats.var_mean_max
    @test size(ref_stats.var_std_max) == (n_cases, n_variables)
    @test ref_stats.centered == false
    @test ref_stats.case_numbers == [case_number]
    @test ref_stats.n_variables == n_variables
    @test ref_stats.n_heights == n_heights
    @test ref_stats.n_times == n_times

    #setup
    stats_config = Dict(
        "normalization" => "std_normalized",
        "perform_pca" => false,
        "variance_loss" => 0.1,
        "tikhonov_mode" => "absolute",
        "tikhonov_noise" => 0.01,
    )

    #action
    ref_stats = KCP.ReferenceStatistics(obs, stats_config, case_number, n_variables, n_heights, n_times)

    #test
    @test ref_stats.y_norm == ref_stats.var_std_max
    @test ref_stats.centered == true
end

@testset "normalize observations" begin
    # setup
    n_cases = 2
    n_heights = 2
    n_times = 2
    n_vars = 2
    n_samples = 20
    obs = rand(n_cases * n_vars * n_heights * n_times, n_samples) .* 10.0
    obs_m = mean(obs, dims = 2)
    norm_vecs = [
        maximum([obs_m[1:2]; obs_m[5:6]]) maximum([obs_m[3:4]; obs_m[7:8]])
        maximum([obs_m[9:10]; obs_m[13:14]]) maximum([obs_m[11:12]; obs_m[15:16]])
    ]
    norm_vecs_extended = reshape(repeat(norm_vecs, inner = [2 2])', :, 1)
    obs_m_normalized = obs_m ./ norm_vecs_extended
    stats_config = Dict(
        "normalization" => "mean_normalized",
        "perform_pca" => false,
        "variance_loss" => 0.1,
        "tikhonov_mode" => "absolute",
        "tikhonov_noise" => 0.01,
    )

    #action
    ref_stats = KCP.combine_ref_stats(KCP.make_ref_stats_list(obs, stats_config, n_cases, n_vars, n_heights, n_times))

    #test
    @test ref_stats.y_norm == ref_stats.var_mean_max == norm_vecs
    @test ref_stats.y_full ≈ obs_m_normalized rtol = eps(Float64) * 10.0
end
