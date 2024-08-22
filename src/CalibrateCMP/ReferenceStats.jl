
Base.@kwdef struct ReferenceStatistics{FT <: Real}
    y::Vector{FT}
    Γ::AbstractMatrix{FT}
    P_pca::Union{AbstractMatrix{FT}, UniformScaling}
    y_full::Vector{FT}
    Γ_full::AbstractMatrix{FT}
    obs_mean::Vector{FT}
    y_norm::Matrix{FT}
    var_mean_max::Matrix{FT}
    var_std_max::Matrix{FT}
    centered::Bool
    case_numbers::Vector{Int}
    n_heights::Vector{Int}
    n_times::Vector{Int}


    function ReferenceStatistics(
        y::Vector{FT},
        Γ::AbstractMatrix{FT},
        P_pca::Union{AbstractMatrix{FT}, UniformScaling},
        y_full::Vector{FT},
        Γ_full::AbstractMatrix{FT},
        obs_mean::Vector{FT},
        y_norm::Matrix{FT},
        var_mean_max::Matrix{FT},
        var_std_max::Matrix{FT},
        centered::Bool,
        case_numbers::Vector{Int},
        n_heights::Vector{Int},
        n_times::Vector{Int},
    )
        return new{FT}(
            y,
            Γ,
            P_pca,
            y_full,
            Γ_full,
            obs_mean,
            y_norm,
            var_mean_max,
            var_std_max,
            centered,
            case_numbers,
            n_heights,
            n_times,
        )
    end

    # single case ref stat
    function ReferenceStatistics(
        obs::Matrix{FT},
        stats_config::Dict,
        case_number::Int,
        n_heights::Vector{Int},
        n_times::Int,
    ) where {FT <: Real}

        @assert sum(n_heights) * n_times == size(obs)[1]
        case_numbers = [case_number]
        _n_cases = 1
        (var_mean_max, var_std_max) = find_mean_and_std_maximum(obs, _n_cases, n_heights, [n_times])
        if stats_config["normalization"] == "mean_normalized"
            y_norm = var_mean_max
            centered = false
        elseif stats_config["normalization"] == "std_normalized"
            y_norm = var_std_max
            centered = true
        else
            error("normalization not recognized!!!")
        end
        _n_variables = length(n_heights)
        weights =
            "weights" in keys(stats_config) ? reshape(stats_config["weights"], 1, _n_variables) : ones(1, _n_variables)
        y_norm = y_norm ./ weights
        obs_normalized = normalize_obs(obs, y_norm, _n_cases, n_heights, [n_times], centered = centered)
        obs_mean = mean(obs, dims = 2)[:]
        y_full = mean(obs_normalized, dims = 2)[:]
        Γ_full = cov(obs_normalized, dims = 2)

        if stats_config["perform_pca"]
            y, Γ, P_pca = obs_PCA(y_full, Γ_full, stats_config["variance_loss"])
        else
            y = y_full
            Γ = Γ_full
            P_pca = 1.0I(length(y))
        end

        # Condition global covariance matrix, PCA
        if stats_config["tikhonov_mode"] == "relative"
            tikhonov_noise = max(stats_config["tikhonov_noise"], 10 * sqrt(eps(FT)))
            Γ = Γ + tikhonov_noise * eigmax(Γ) * I
        else
            Γ = Γ + stats_config["tikhonov_noise"] * I
        end

        isposdef(Γ) ? nothing : println("WARNING: Covariance matrix Γ is ill-conditioned, consider regularization.")

        return new{FT}(
            y,
            Γ,
            P_pca,
            y_full,
            Γ_full,
            obs_mean,
            y_norm,
            var_mean_max,
            var_std_max,
            centered,
            case_numbers,
            n_heights,
            [n_times],
        )
    end

end

function find_mean_and_std_maximum(
    y::Matrix{FT},
    n_cases::Int,
    n_heights::Vector{Int},
    n_times::Vector{Int},
) where {FT <: Real}

    _n_variables = length(n_heights)
    _n_single_case_sim = sum(n_heights) * n_times
    @assert size(y)[1] == sum(_n_single_case_sim)

    var_mean_max = ones(n_cases, _n_variables)
    var_std_max = ones(n_cases, _n_variables)
    for i in 1:n_cases
        _y_single_case = y[(sum(_n_single_case_sim[1:(i - 1)]) + 1):sum(_n_single_case_sim[1:i]), :]

        _y_mean_vector_single_case = mean(_y_single_case, dims = 2)
        _y_mean_matrix::Matrix{FT} = reshape(_y_mean_vector_single_case, sum(n_heights), n_times[i])
        _var_mean_max = [
            maximum(abs.(_y_mean_matrix[(sum(n_heights[1:(j - 1)]) + 1):sum(n_heights[1:j]), :])) for
            j in 1:_n_variables
        ]
        _var_mean_max = ifelse.(_var_mean_max .> eps(FT), _var_mean_max, 1.0)
        var_mean_max[i, :] = _var_mean_max

        _y_std_vector_single_case = std(_y_single_case, dims = 2)
        _y_std_matrix::Matrix{FT} = reshape(_y_std_vector_single_case, sum(n_heights), n_times[i])
        _var_std_max = [
            maximum(abs.(_y_std_matrix[(sum(n_heights[1:(j - 1)]) + 1):sum(n_heights[1:j]), :])) for j in 1:_n_variables
        ]
        _var_std_max = ifelse.(_var_std_max .> eps(FT), _var_std_max, 1.0)
        var_std_max[i, :] = _var_std_max

    end
    return (var_mean_max = var_mean_max, var_std_max = var_std_max)
end

function normalize_obs(
    y::Matrix{FT},
    ynorm::Matrix{FT},
    n_cases::Int,
    n_heights::Vector{Int},
    n_times::Vector{Int};
    centered::Bool = true,
) where {FT <: Real}

    _n_variables = length(n_heights)
    _n_single_case_sim = sum(n_heights) * n_times
    @assert size(y)[1] == sum(_n_single_case_sim)
    @assert size(ynorm) == (n_cases, _n_variables)

    _y_normalized = centered ? y .- mean(y, dims = 2) : copy(y)
    for i in 1:n_cases
        # normalise part of y belonging to case i by _norm_vec
        _norm_vec = ynorm[i, :]
        _norm_vec_extended_single_time = zeros(sum(n_heights))
        for j in 1:_n_variables
            _norm_vec_extended_single_time[(sum(n_heights[1:(j - 1)]) + 1):sum(n_heights[1:j])] .= _norm_vec[j]
        end
        _norm_vec_extended = reshape(ones(sum(n_heights), n_times[i]) .* _norm_vec_extended_single_time, :, 1)
        _row_indeces = (sum(_n_single_case_sim[1:(i - 1)]) + 1):sum(_n_single_case_sim[1:i])
        _y_normalized[_row_indeces, :] = _y_normalized[_row_indeces, :] ./ _norm_vec_extended
    end

    return _y_normalized
end

function normalize_sim(
    sim_vec::Vector{FT},
    ref_stats::ReferenceStatistics;
    norm::Union{Matrix{FT}, Nothing} = nothing,
) where {FT <: Real}

    _n_cases = length(ref_stats.case_numbers)
    _n_heights = ref_stats.n_heights
    _n_variables = length(_n_heights)
    _n_times = ref_stats.n_times
    _ynorm = norm === nothing ? ref_stats.y_norm : norm

    _sim_vec_normalized = ref_stats.centered ? sim_vec .- ref_stats.obs_mean : copy(sim_vec)
    _n_single_case_sim = sum(_n_heights) * _n_times
    for i in 1:_n_cases
        # normalise part of y belonging to case i by _norm_vec
        _norm_vec = _ynorm[i, :]
        _norm_vec_extended_single_time = zeros(sum(_n_heights))
        for j in 1:_n_variables
            _norm_vec_extended_single_time[(sum(_n_heights[1:(j - 1)]) + 1):sum(_n_heights[1:j])] .= _norm_vec[j]
        end
        _norm_vec_extended = reshape(ones(sum(_n_heights), _n_times[i]) .* _norm_vec_extended_single_time, :, 1)
        _indeces = (sum(_n_single_case_sim[1:(i - 1)]) + 1):sum(_n_single_case_sim[1:i])
        _sim_vec_normalized[_indeces] = _sim_vec_normalized[_indeces] ./ _norm_vec_extended
    end

    return _sim_vec_normalized
end

function unnormalize_sim(
    sim_vec::Vector{FT},
    ref_stats::ReferenceStatistics;
    norm::Union{Matrix{FT}, Nothing} = nothing,
) where {FT <: Real}

    _n_cases = length(ref_stats.case_numbers)
    _n_heights = ref_stats.n_heights
    _n_variables = length(_n_heights)
    _n_times = ref_stats.n_times
    _ynorm = norm === nothing ? ref_stats.y_norm : norm

    _sim_vec_unnormalized = copy(sim_vec)
    _n_single_case_sim = sum(_n_heights) * _n_times
    for i in 1:_n_cases
        # unnormalise part of y belonging to case i
        _norm_vec = _ynorm[i, :]
        _norm_vec_extended_single_time = zeros(sum(_n_heights))
        for j in 1:_n_variables
            _norm_vec_extended_single_time[(sum(_n_heights[1:(j - 1)]) + 1):sum(_n_heights[1:j])] .= _norm_vec[j]
        end
        _norm_vec_extended = reshape(ones(sum(_n_heights), _n_times[i]) .* _norm_vec_extended_single_time, :, 1)
        _indeces = (sum(_n_single_case_sim[1:(i - 1)]) + 1):sum(_n_single_case_sim[1:i])
        _sim_vec_unnormalized[_indeces] = _sim_vec_unnormalized[_indeces] .* _norm_vec_extended
    end
    if ref_stats.centered
        _sim_vec_unnormalized = _sim_vec_unnormalized .+ ref_stats.obs_mean
    end

    return _sim_vec_unnormalized
end

function obs_PCA(y_mean::Vector{FT}, y_var::AbstractMatrix{FT}, allowed_var_loss::FT = 1.0e-1) where {FT <: Real}
    λ_pca, P_pca = pca(y_var, allowed_var_loss)
    # Project mean
    y_pca = P_pca' * y_mean
    y_var_pca = Diagonal(λ_pca)
    return y_pca, y_var_pca, P_pca
end

function pca(covmat::AbstractMatrix{FT}, allowed_var_loss::FT) where {FT <: Real}
    eigvals, eigvecs = eigen(covmat)
    # Get index of leading eigenvalues, eigvals are ordered from low to high in julia
    # This expression recovers 1 extra eigenvalue compared to threshold
    leading_eigs = findall(<(1.0 - allowed_var_loss), -cumsum(eigvals) / sum(eigvals) .+ 1)
    P_pca = eigvecs[:, leading_eigs]
    λ_pca = eigvals[leading_eigs]
    return λ_pca, P_pca
end

function pca_transform(vec::Vector{FT}, ref_stats::ReferenceStatistics) where {FT <: Real}
    return ref_stats.P_pca' * vec
end

function pca_itransform(vec::Vector{FT}, ref_stats::ReferenceStatistics) where {FT <: Real}
    return ref_stats.P_pca * vec
end

# returns a list of ref stats for all cases, each for a single case
function make_ref_stats_list(
    obs::Matrix{FT},
    stats_config::Dict,
    n_cases::Int,
    n_heights::Vector{Int},
    n_times::Vector{Int},
) where {FT <: Real}

    _n_single_case_sim = sum(n_heights) * n_times
    @assert sum(_n_single_case_sim) == size(obs)[1]
    ref_stats_list = ReferenceStatistics[]
    for i in 1:n_cases
        _single_case_row_indices = (sum(_n_single_case_sim[1:(i - 1)]) + 1):sum(_n_single_case_sim[1:i])
        _obs_single_case = obs[_single_case_row_indices, :]

        ref_stats_list = [
            ref_stats_list
            ReferenceStatistics(_obs_single_case, stats_config, i, n_heights, n_times[i])
        ]
    end
    return ref_stats_list
end

function combine_ref_stats(ref_stats_list::Vector{ReferenceStatistics})
    @assert !isempty(ref_stats_list)

    y = ref_stats_list[1].y
    Γ = ref_stats_list[1].Γ
    P_pca = ref_stats_list[1].P_pca
    y_full = ref_stats_list[1].y_full
    Γ_full = ref_stats_list[1].Γ_full
    obs_mean = ref_stats_list[1].obs_mean
    y_norm = ref_stats_list[1].y_norm
    var_mean_max = ref_stats_list[1].var_mean_max
    var_std_max = ref_stats_list[1].var_std_max
    centered = ref_stats_list[1].centered
    case_numbers = ref_stats_list[1].case_numbers
    n_heights = ref_stats_list[1].n_heights
    n_times = ref_stats_list[1].n_times

    for ref_stats in ref_stats_list[2:end]
        @assert centered == ref_stats.centered
        @assert n_heights == ref_stats.n_heights
        n_times = [n_times; ref_stats.n_times]
        case_numbers = [case_numbers; ref_stats.case_numbers]
        y_norm = [y_norm; ref_stats.y_norm]
        var_mean_max = [var_mean_max; ref_stats.var_mean_max]
        var_std_max = [var_std_max; ref_stats.var_std_max]
        obs_mean = [obs_mean; ref_stats.obs_mean]
        y_full = [y_full; ref_stats.y_full]
        y = [y; ref_stats.y]
        P_pca = make_block_diagonal_matrix(P_pca, ref_stats.P_pca)
        Γ_full = make_block_diagonal_matrix(Γ_full, ref_stats.Γ_full)
        Γ = make_block_diagonal_matrix(Γ, ref_stats.Γ)
    end

    RS = ReferenceStatistics(;
        y = y,
        Γ = Γ,
        P_pca = P_pca,
        y_full = y_full,
        Γ_full = Γ_full,
        obs_mean = obs_mean,
        y_norm = y_norm,
        var_mean_max = var_mean_max,
        var_std_max = var_std_max,
        centered = centered,
        case_numbers = case_numbers,
        n_heights = n_heights,
        n_times = n_times,
    )
    return RS
end
