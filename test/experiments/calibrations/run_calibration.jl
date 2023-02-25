import Kinematic1D
const KiD = Kinematic1D

include("./config.jl")

data_save_directory = KiD.make_output_directories()

config = get_config()

priors = KiD.construct_priors(config["prior"]["parameters"])

obs = KiD.get_obs!(config)
ref_stats_list = KiD.make_ref_stats_list(obs, config["statistics"], KiD.get_numbers_from_config(config)...)

u_names = keys(config["prior"]["parameters"])

if config["process"]["method"] == "EKP"
    res = KiD.calibrate(KiD.EKPStyle(), priors, config, ref_stats_list)
    KiD.ensemble_convergence(res, priors, config, file_name = data_save_directory * "ensemble_convergence.gif")
elseif config["process"]["method"] == "Optim"
    res = KiD.calibrate(KiD.OptimStyle(), priors, config, KiD.combine_ref_stats(ref_stats_list))
elseif config["process"]["method"] == "both"
    res = KiD.calibrate(KiD.EKPStyle(), priors, config, ref_stats_list)
    u_bests = KiD.get_results(res, priors).ϕ_optim
    for (i, k) in enumerate(u_names)
        v = config["prior"]["parameters"][k]
        config["prior"]["parameters"][k] = (mean = u_bests[i], var = v.var, lbound = v.lbound, ubound = v.ubound)
    end
    res = KiD.calibrate(KiD.OptimStyle(), priors, config, KiD.combine_ref_stats(ref_stats_list))
else
    error("Calibration method not implemented!")
end

ϕ_bests = KiD.get_results(res, priors)
KiD.print_results(ϕ_bests.ϕ_optim, u_names)
KiD.save_data(
    res,
    ϕ_bests,
    collect(u_names),
    config,
    file_name = data_save_directory * config["process"]["output_file_name"],
)
