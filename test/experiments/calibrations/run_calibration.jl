import KinematicDriver.CalibrateCMP as KCP

include("./config.jl")

data_save_directory = KCP.make_output_directories(dir = joinpath(@__DIR__, "/run_calibration_output/"))

config = get_config()

priors = KCP.construct_priors(config["prior"]["parameters"])

obs = KCP.get_obs!(config)
ref_stats_list = KCP.make_ref_stats_list(obs, config["statistics"], KCP.get_numbers_from_config(config)...)

u_names = keys(config["prior"]["parameters"])

if config["process"]["method"] == "EKP"
    res = KCP.calibrate(KCP.EKPStyle(), priors, config, ref_stats_list)
    KCP.ensemble_convergence(res, priors, config, file_name = data_save_directory * "ensemble_convergence.gif")
elseif config["process"]["method"] == "Optim"
    res = KCP.calibrate(KCP.OptimStyle(), priors, config, KCP.combine_ref_stats(ref_stats_list))
elseif config["process"]["method"] == "both"
    res = KCP.calibrate(KCP.EKPStyle(), priors, config, ref_stats_list)
    u_bests = KCP.get_results(res, priors).ϕ_optim
    for (i, k) in enumerate(u_names)
        v = config["prior"]["parameters"][k]
        config["prior"]["parameters"][k] = (mean = u_bests[i], var = v.var, lbound = v.lbound, ubound = v.ubound)
    end
    res = KCP.calibrate(KCP.OptimStyle(), priors, config, KCP.combine_ref_stats(ref_stats_list))
else
    error("Calibration method not implemented!")
end

ϕ_bests = KCP.get_results(res, priors)
KCP.print_results(ϕ_bests.ϕ_optim, u_names)
KCP.save_data(
    res,
    ϕ_bests,
    collect(u_names),
    config,
    file_name = data_save_directory * config["process"]["output_file_name"],
)
