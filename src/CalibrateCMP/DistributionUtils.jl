function initial_parameter_ensemble(priors, N_ensemble; rng_seed = 10)
    Random.seed!(rng_seed)
    initial_ensemble = construct_initial_ensemble(priors, N_ensemble)
    return initial_ensemble
end

function construct_priors(
    params::Dict{String, @NamedTuple{mean::FT, var::FT, lbound::FT, ubound::FT}},
) where {FT <: Real}
    u_names = collect(keys(params))
    u_values = collect(values(params))
    n_params = length(u_names)

    # All vars are approximated as Gaussian in unconstrained space
    marginal_priors = Vector{ParameterDistribution}(undef, n_params)
    for i in 1:n_params
        marginal_priors[i] =
            constrained_gaussian(u_names[i], u_values[i].mean, u_values[i].var, u_values[i].lbound, u_values[i].ubound)
    end
    return combine_distributions(marginal_priors)
end
