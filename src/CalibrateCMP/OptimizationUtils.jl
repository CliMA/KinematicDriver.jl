abstract type AbstractOptimizationStyle end
struct EKPStyle <: AbstractOptimizationStyle end
struct OptimStyle <: AbstractOptimizationStyle end

function compute_loss(u::Array{Float64, 1}, u_names::Array{String, 1}, config::Dict, RS::ReferenceStatistics)
    G_n = run_dyn_model(u, u_names, config, RS = RS)
    yg_diff = G_n - RS.y
    n_cases = length(RS.case_numbers)
    loss = dot(yg_diff, RS.Γ \ yg_diff) / n_cases
    return loss
end

function calibrate(os::AbstractOptimizationStyle, priors, config, RS::ReferenceStatistics; verbose::Int)
    error("calibrate not implemented for a given $os")
end

# Optim
function calibrate(::OptimStyle, priors, config, RS::ReferenceStatistics; verbose::Int = 1)
    _u_names = collect(keys(config["prior"]["parameters"]))
    _u_inits = collect([v.mean for v in values(config["prior"]["parameters"])])
    _θ_inits = mapslices(x -> transform_constrained_to_unconstrained(priors, x), _u_inits; dims = 1)
    _f = make_optim_loss_function(_u_names, priors, config, RS)

    return optimize(
        _f,
        _θ_inits;
        g_abstol = config["process"]["tol"],
        iterations = config["process"]["maxiter"],
        show_trace = verbose > 0,
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
function calibrate(::EKPStyle, priors, config, ref_stats_list::Vector{ReferenceStatistics}; verbose::Int = 3)
    @assert config["process"]["EKP_method"] == "EKI" || config["process"]["EKP_method"] == "UKI" || config["process"]["EKP_method"] == "ETKI"

    u_names = collect(keys(config["prior"]["parameters"]))
    n_iter = config["process"]["n_iter"]
    if config["process"]["EKP_method"] == "EKI" || config["process"]["EKP_method"] == "ETKI"
        n_ensemble = config["process"]["n_ens"]
    else
        n_ensemble = 2 * length(u_names) + 1
    end
    param_ensemble = [initial_parameter_ensemble(priors, n_ensemble)]
    n_cases = length(ref_stats_list)
    batch_size = config["process"]["batch_size"]
    n_batches = floor(Int, n_cases / batch_size)
    for n in 1:n_iter

        if verbose >= 1
            println("==========================================================================================")
            println("iteration: ", n)
            if verbose >= 3
                RS = combine_ref_stats(ref_stats_list)
                ϕ_mean = get_ϕ_mean_final(param_ensemble, priors)
                model_error = compute_error_metrics(ϕ_mean, u_names, config, RS)
                println(
                    "loss = ",
                    model_error.loss,
                    ",\t mse_m = ",
                    model_error.mse_m,
                    ",\t mse_s = ",
                    model_error.mse_s,
                )
            end
        end

        cases_indices_shuffled = shuffle!(collect(1:n_cases))
        for j in 1:n_batches
            batch_indices = cases_indices_shuffled[((j - 1) * batch_size + 1):(j * batch_size)]
            RS = combine_ref_stats(ref_stats_list[batch_indices])

            ϕ_mean = get_ϕ_mean_final(param_ensemble, priors)
            if verbose >= 1
                println("---------------------------")
                println("batch: ", j, ", cases: ", batch_indices)
                for (i, name) in enumerate(u_names)
                    print(name, " = ", round(ϕ_mean[i], sigdigits = 4), ", ")
                end
                print("\n")
                if verbose >= 2
                    model_error = compute_error_metrics(ϕ_mean, u_names, config, RS)
                    println(
                        "loss = ",
                        model_error.loss,
                        ",\t mse_m = ",
                        model_error.mse_m,
                        ",\t mse_s = ",
                        model_error.mse_s,
                    )
                end
            end

            ekpobj = generate_ekp(param_ensemble[end], RS, priors, config["process"])
            ϕ_n = get_ϕ_ensemble(param_ensemble[end], priors)
            G_n = run_ensembles(ϕ_n, u_names, config, n_ensemble, RS = RS)
            G_ens = hcat(G_n...)
            if config["process"]["augmented"]
                G_ens = vcat(G_ens, param_ensemble[end])
            end
            EnsembleKalmanProcesses.update_ensemble!(ekpobj, G_ens)
            param_ensemble = [param_ensemble; [get_u_final(ekpobj)]]
        end
    end

    return param_ensemble
end

function generate_ekp(
    param_ensemble::Matrix{FT},
    RS::ReferenceStatistics,
    priors::ParameterDistribution,
    process_settings,
) where {FT <: Real}

    if process_settings["augmented"]
        μ = vcat(mean(priors)...)
        Γ_θ = cov(priors)
        y = vcat(RS.y, μ)
        Γ = cat(RS.Γ, Γ_θ .* process_settings["reg_degree_scale"], dims = (1, 2))
    else
        y = RS.y
        Γ = RS.Γ
    end

    if process_settings["EKP_method"] == "EKI"
        ekpobj = EnsembleKalmanProcess(
            param_ensemble,
            y,
            Γ,
            Inversion(),
            scheduler = DefaultScheduler(process_settings["Δt"]),
        )
    elseif process_settings["EKP_method"] == "ETKI"
        ekpobj = EnsembleKalmanProcess(
            param_ensemble,
            y,
            Γ,
            TransformInversion(),
            scheduler = DefaultScheduler(process_settings["Δt"]),
        )
    elseif process_settings["EKP_method"] == "UKI"
        process = Unscented(
            mean(param_ensemble, dims = 2)[:],
            cov(param_ensemble, dims = 2);
            α_reg = process_settings["α_reg"],
            update_freq = process_settings["update_freq"],
        )
        ekpobj = EnsembleKalmanProcess(y, Γ, process, scheduler = DefaultScheduler(process_settings["Δt"]))
    end
    return ekpobj
end

function get_ϕ_ensemble(param_ensemble::Matrix{FT}, priors::ParameterDistribution) where {FT <: Real}
    _ϕ_ensemble = zeros(size(param_ensemble))
    for j in 1:size(_ϕ_ensemble)[2]
        _ϕ_ensemble[:, j] = transform_unconstrained_to_constrained(priors, param_ensemble[:, j])
    end
    return _ϕ_ensemble
end

function get_ϕ_mean_final(param_ensemble::Vector{Matrix{FT}}, priors::ParameterDistribution) where {FT <: Real}
    _ϕ_mean = transform_unconstrained_to_constrained(priors, mean(param_ensemble[end], dims = 2))[:]
    return _ϕ_mean
end

function get_results(param_ensemble::Vector{Matrix{FT}}, priors::ParameterDistribution) where {FT <: Real}
    _ϕ_optim = get_ϕ_mean_final(param_ensemble, priors)
    _cov_optim = cov(param_ensemble[end], dims = 2)
    _cov_optim = _cov_optim ./ tr(_cov_optim) .* length(_ϕ_optim)
    return (ϕ_optim = _ϕ_optim, cov_optim = _cov_optim)
end

function run_ensembles(
    u::Array{FT, 2},
    u_names::Array{String, 1},
    config::Dict,
    n_ensemble::Int;
    RS::Union{Nothing, ReferenceStatistics} = nothing,
    case_numbers::Vector{Int} = Int[],
)
    G_n = zeros(length(RS.y), n_ensemble)
    Threads.@threads for i in 1:n_ensemble
        # run the model with the current parameters, i.e., map θ to G(θ)
        G_n[:, i] = run_dyn_model(u[:, i], u_names, config, RS = RS, case_numbers=case_numbers) 
    end
    return [G_n[:,i] for i in 1:n_ensemble]
end
