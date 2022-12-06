abstract type AbstractOptimizationStyle end
struct EKPStyle <: AbstractOptimizationStyle end
struct OptimStyle <: AbstractOptimizationStyle end

function compute_loss(u::Array{Float64, 1}, u_names::Array{String, 1}, config::Dict, RS::ReferenceStatistics)
    G_t = run_dyn_model(u, u_names, config, RS = RS)
    yg_diff_pca = G_t - RS.y
    loss = dot(yg_diff_pca, RS.Γ \ yg_diff_pca)
    return loss
end

function calibrate(os::AbstractOptimizationStyle, priors, config, RS::ReferenceStatistics; verbose::Bool)
    error("calibrate not implemented for a given $os")
end

# Optim
function calibrate(::OptimStyle, priors, config, RS::ReferenceStatistics; verbose::Bool = true)
    _u_names = collect(keys(config["prior"]["parameters"]))
    _u_inits = collect([v.mean for v in values(config["prior"]["parameters"])])
    _θ_inits = mapslices(x -> transform_constrained_to_unconstrained(priors, x), _u_inits; dims = 1)
    _f = make_optim_loss_function(_u_names, priors, config, RS)

    return optimize(
        _f,
        _θ_inits;
        g_abstol = config["process"]["tol"],
        iterations = config["process"]["maxiter"],
        show_trace = verbose,
    )
end

function make_optim_loss_function(u_names::Array{String, 1}, priors, config::Dict, RS::ReferenceStatistics)
    function optim_loss_function(θ::Array{Float64, 1})
        u = mapslices(x -> transform_unconstrained_to_constrained(priors, x), θ; dims = 1)
        return compute_loss(u, u_names, config, RS)
    end
    return optim_loss_function
end

function get_results(optres::Optim.MultivariateOptimizationResults, priors::ParameterDistribution)
    _θ_optim = Optim.minimizer(optres)
    _ϕ_optim = mapslices(x -> transform_unconstrained_to_constrained(priors, x), _θ_optim; dims = 1)
    _cov_optim = Matrix(I(length(_ϕ_optim))) .* eltype(_ϕ_optim)(1)
    return (ϕ_optim = _ϕ_optim, cov_optim = _cov_optim)
end

# EKP
function calibrate(::EKPStyle, priors, config, RS::ReferenceStatistics; verbose::Bool = true)
    method = config["process"]["EKP_method"]
    if method == "EKI"
        ekpobj = EnsembleKalmanProcess(
            initial_parameter_ensemble(priors, config["process"]["n_ens"]),
            RS.y,
            RS.Γ,
            Inversion(),
            Δt = config["process"]["Δt"],
        )
    elseif method == "UKI"
        α_reg = config["process"]["α_reg"]
        update_freq = config["process"]["update_freq"]
        prior_mean = EKP.ParameterDistributions.mean(priors)
        prior_cov = EKP.ParameterDistributions.cov(priors)
        process = Unscented(prior_mean, prior_cov; α_reg = α_reg, update_freq = update_freq)
        ekpobj = EnsembleKalmanProcess(RS.y, RS.Γ, process, Δt = config["process"]["Δt"])
    else
        error("EKP method not implemented!!!")
    end


    u_names = collect(keys(config["prior"]["parameters"]))
    n_iter_min = config["process"]["n_iter_min"]
    n_iter_max = config["process"]["n_iter_max"]
    n_ensemble = config["process"]["EKP_method"] == "EKI" ? config["process"]["n_ens"] : 2 * length(u_names) + 1
    old_loss = Float64(1e16)
    for n in 1:n_iter_max
        ϕ_n = get_ϕ_final(priors, ekpobj)
        ϕ_mean = mean(ϕ_n, dims = 2)
        loss = compute_loss(ϕ_mean[:], u_names, config, RS)
        if verbose
            println("iteration ", n, " : loss = ", loss)
            println(round.(ϕ_mean, sigdigits = 4))
        end
        if (loss > old_loss) & (n > n_iter_min)
            break
        end
        old_loss = loss

        G_n = [run_dyn_model(ϕ_n[:, i], u_names, config, RS = RS) for i in 1:n_ensemble]
        G_ens = hcat(G_n...)
        EnsembleKalmanProcesses.update_ensemble!(ekpobj, G_ens)
    end

    return ekpobj
end

function get_results(ekpobj::EnsembleKalmanProcess, priors::ParameterDistribution)
    _ϕ_optim = get_ϕ_mean_final(priors, ekpobj)
    _cov_optim = get_u_cov_final(ekpobj)
    _cov_optim = _cov_optim ./ tr(_cov_optim) .* length(_ϕ_optim)
    return (ϕ_optim = _ϕ_optim, cov_optim = _cov_optim)
end
