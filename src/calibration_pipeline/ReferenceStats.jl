
Base.@kwdef struct ReferenceStatistics{FT <: Real}
    y::Vector{FT}
    Γ::AbstractMatrix{FT}
    P_pca::Union{AbstractMatrix{FT}, UniformScaling}
    y_full::Vector{FT}
    Γ_full::AbstractMatrix{FT}

    function ReferenceStatistics(truth::Observation{FT}, stats_config::Dict) where {FT <: Real}

        y_full, Γ_full = get_obs(truth)
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

        return new{FT}(y, Γ, P_pca, y_full, Γ_full)
    end

end

function get_obs(truth::Observation{FT}) where {FT <: Real}
    # Get true observables
    y_ = truth.mean
    Γ_ = truth.obs_noise_cov
    return y_, Γ_
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
