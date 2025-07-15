TT.@testset "Constructing priors" begin
    #setup
    priors_config = get_prior_config()
    param_names = keys(priors_config["parameters"])

    #action
    priors = construct_priors(priors_config["parameters"])

    #test
    TT.@test priors isa EKP.ParameterDistribution{EKP.Parameterized, EKP.ParameterDistributions.Constraint, String}
    TT.@test Set(priors.name) == Set(param_names)
    TT.@test length(priors.constraint) == length(priors.distribution) == length(param_names)
end

TT.@testset "Initial parameter ensemble" begin
    #setup
    priors_config = get_prior_config()
    priors = construct_priors(priors_config["parameters"])
    n_params = length(keys(priors_config["parameters"]))
    n_ensemble = 20

    #action
    ensemble = initial_parameter_ensemble(priors, n_ensemble, rng_seed = 15)

    #test
    TT.@test ensemble isa Matrix
    TT.@test size(ensemble) == (n_params, n_ensemble)
end
