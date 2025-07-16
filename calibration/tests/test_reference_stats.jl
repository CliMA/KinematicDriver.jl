TT.@testset "Test PCA computations" begin
    #setup
    y_mean = rand(10)
    covmat = convert(Array, LA.Diagonal(Random.shuffle(collect(range(0.1, 1.0, 10)))))
    allowed_var_loss = 0.1

    #action
    λ_pca, P_pca = pca(covmat, allowed_var_loss)

    #test
    TT.@test λ_pca == collect(range(0.3, 1.0, 8))
    TT.@test maximum(P_pca, dims = 1) == maximum(P_pca, dims = 1) == ones(1, 8)

    #action
    y_pca, covmat_pca, P_pca = obs_PCA(y_mean, covmat, allowed_var_loss)

    #action
    TT.@test length(y_pca) == 8
    TT.@test y_pca == P_pca' * y_mean
    TT.@test covmat_pca == LA.Diagonal(λ_pca)
    TT.@test maximum(P_pca, dims = 1) == maximum(P_pca, dims = 1) == ones(1, 8)
end

TT.@testset "Constructing reference statistics" begin
    #setup
    stats_config = Dict(
        "normalization" => "mean_normalized",
        "perform_pca" => true,
        "variance_loss" => 0.1,
        "tikhonov_mode" => "absolute",
        "tikhonov_noise" => 0.0,
    )
    n_cases = 1
    n_heights = [3, 3]
    n_variables = length(n_heights)
    n_times = [2]
    n_tot = sum(sum(n_heights) * n_times)
    n_samples = 50
    case_number = 1
    obs = rand(n_tot, n_samples)
    covmat = Statistics.cov(obs, dims = 2)

    #action
    ref_stats = ReferenceStatistics(obs, stats_config, case_number, n_heights, n_times[1])

    #test
    TT.@test length(ref_stats.y_full) == n_tot
    TT.@test size(ref_stats.Γ_full) == (n_tot, n_tot)
    TT.@test length(ref_stats.y) < n_tot
    TT.@test sum(ref_stats.Γ, dims = 1)[:] == LA.diag(ref_stats.Γ)
    TT.@test size(ref_stats.Γ)[1] == size(ref_stats.Γ)[2] < n_tot
    TT.@test size(ref_stats.P_pca)[1] == n_tot
    TT.@test size(ref_stats.P_pca)[2] < n_tot
    TT.@test ref_stats.obs_mean ≈ Statistics.mean(obs, dims = 2) rtol = eps(Float64) * 10.0
    TT.@test size(ref_stats.y_norm) == (n_cases, n_variables)
    TT.@test ref_stats.y_norm == ref_stats.var_mean_max
    TT.@test size(ref_stats.var_std_max) == (n_cases, n_variables)
    TT.@test ref_stats.centered == false
    TT.@test ref_stats.case_numbers == [case_number]
    TT.@test ref_stats.n_heights == n_heights
    TT.@test ref_stats.n_times == n_times

    #setup
    stats_config = Dict(
        "normalization" => "mean_normalized",
        "perform_pca" => false,
        "variance_loss" => 0.1,
        "tikhonov_mode" => "absolute",
        "tikhonov_noise" => 0.01,
    )

    #action
    ref_stats = ReferenceStatistics(obs, stats_config, case_number, n_heights, n_times[1])

    #test
    TT.@test length(ref_stats.y_full) == n_tot
    TT.@test size(ref_stats.Γ_full) == (n_tot, n_tot)
    TT.@test ref_stats.y == ref_stats.y_full
    TT.@test ref_stats.Γ == ref_stats.Γ_full .+ 0.01 * LA.I(n_tot)
    TT.@test ref_stats.P_pca == LA.I(n_tot)
    TT.@test ref_stats.obs_mean ≈ Statistics.mean(obs, dims = 2) rtol = eps(Float64) * 10.0
    TT.@test size(ref_stats.y_norm) == (n_cases, n_variables)
    TT.@test ref_stats.y_norm == ref_stats.var_mean_max
    TT.@test size(ref_stats.var_std_max) == (n_cases, n_variables)
    TT.@test ref_stats.centered == false
    TT.@test ref_stats.case_numbers == [case_number]
    TT.@test ref_stats.n_heights == n_heights
    TT.@test ref_stats.n_times == n_times

    #setup
    stats_config = Dict(
        "normalization" => "std_normalized",
        "perform_pca" => false,
        "variance_loss" => 0.1,
        "tikhonov_mode" => "absolute",
        "tikhonov_noise" => 0.01,
    )

    #action
    ref_stats = ReferenceStatistics(obs, stats_config, case_number, n_heights, n_times[1])

    #test
    TT.@test ref_stats.y_norm == ref_stats.var_std_max
    TT.@test ref_stats.centered == true
end

TT.@testset "normalize observations" begin
    # setup
    n_cases = 2
    n_heights = [2, 2]
    n_times = [2, 2]
    n_samples = 20
    obs = rand(sum(sum(n_heights) * n_times), n_samples) .* 10.0
    obs_m = Statistics.mean(obs, dims = 2)
    norm_vecs = [
        maximum([obs_m[1:2]; obs_m[5:6]]) maximum([obs_m[3:4]; obs_m[7:8]]);
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
    ref_stats = combine_ref_stats(make_ref_stats_list(obs, stats_config, n_cases, n_heights, n_times))

    #test
    TT.@test ref_stats.y_norm == ref_stats.var_mean_max == norm_vecs
    TT.@test ref_stats.y_full ≈ obs_m_normalized rtol = eps(Float64) * 10.0
end
