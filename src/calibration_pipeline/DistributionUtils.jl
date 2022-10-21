function initial_parameter_ensemble(priors, N_ensemble; rng_seed = 10)
    Random.seed!(rng_seed)
    initial_ensemble = construct_initial_ensemble(priors, N_ensemble; rng_seed = rng_seed)
    return initial_ensemble
end

function construct_priors(params::Dict)
    u_names = collect(keys(params))
    u_values = collect(values(params))
    u_means = collect([v.mean for v in u_values])
    u_vars = collect([v.var for v in u_values])
    u_lbounds = collect([v.lbound for v in u_values])
    u_ubounds = collect([v.ubound for v in u_values])

    # All vars are approximated as Gaussian in unconstrained space
    marginal_priors = constrained_gaussian.(u_names, u_means, u_vars, u_lbounds, u_ubounds)
    return combine_distributions(marginal_priors)
end
