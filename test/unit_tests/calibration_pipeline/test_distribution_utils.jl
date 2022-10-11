@testset "Constructing priors" begin
    #setup
    priors_config = get_prior_config()
    param_names = keys(priors_config["parameters"])

    #action
    priors = KID.construct_priors(priors_config["parameters"])

    #test
    @test priors isa ParameterDistribution{Parameterized, Constraint, String}
    @test Set(priors.name) == Set(param_names)
    @test length(priors.constraint) == length(priors.distribution) == length(param_names)
end

@testset "Initial parameter ensemble" begin
    #setup
    priors_config = get_prior_config()
    priors = KID.construct_priors(priors_config["parameters"])
    n_params = length(keys(priors_config["parameters"]))
    n_ensemble = 20

    #action
    ensemble = KID.initial_parameter_ensemble(priors, n_ensemble, rng_seed = 15)

    #test
    @test ensemble isa Matrix
    @test size(ensemble) == (n_params, n_ensemble)
end
