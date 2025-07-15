include("../Calibration.jl")

config_dir = "K1D"
data_save_directory = make_output_directories(joinpath(@__DIR__, config_dir, "calibration_output"))

include(joinpath(@__DIR__, config_dir, "config.jl"))
config = get_config()

priors = construct_priors(config["prior"]["parameters"])

obs = get_obs!(config)
ref_stats_list = make_ref_stats_list(obs, config["statistics"], get_numbers_from_config(config)...)

u_names = keys(config["prior"]["parameters"])

if config["process"]["method"] == "EKP"
    res = calibrate(EKPStyle(), priors, config, ref_stats_list)
    ensemble_convergence(res, priors, config, file_name = joinpath(data_save_directory, "ensemble_convergence.gif"))
elseif config["process"]["method"] == "Optim"
    res = calibrate(OptimStyle(), priors, config, combine_ref_stats(ref_stats_list))
elseif config["process"]["method"] == "both"
    res = calibrate(EKPStyle(), priors, config, ref_stats_list)
    u_bests = get_results(res, priors).ϕ_optim
    for (i, k) in enumerate(u_names)
        v = config["prior"]["parameters"][k]
        config["prior"]["parameters"][k] = (mean = u_bests[i], var = v.var, lbound = v.lbound, ubound = v.ubound)
    end
    res = calibrate(OptimStyle(), priors, config, combine_ref_stats(ref_stats_list))
else
    error("Calibration method not implemented!")
end

ϕ_bests = get_results(res, priors)
print_results(ϕ_bests.ϕ_optim, u_names)
save_data(
    res,
    ϕ_bests,
    collect(u_names),
    config,
    file_name = joinpath(data_save_directory, config["process"]["output_file_name"]),
)
