import Kinematic1D
const KID = Kinematic1D

include("./config.jl")

data_save_directory = KID.make_output_directories()

config = get_config()

priors = KID.construct_priors(config["prior"]["parameters"])

truth = KID.get_obs!(config)
ref_stats = KID.ReferenceStatistics(truth, config["statistics"])

u_names = keys(config["prior"]["parameters"])

if config["process"]["method"] == "EKP"
    res = KID.calibrate(KID.EKPStyle(), priors, config, ref_stats)
    KID.ensemble_convergence(res, priors, config, file_name = data_save_directory * "ensemble_convergence.gif")
elseif config["process"]["method"] == "Optim"
    res = KID.calibrate(KID.OptimStyle(), priors, config, ref_stats)
elseif config["process"]["method"] == "both"
    res = KID.calibrate(KID.EKPStyle(), priors, config, ref_stats)
    u_bests = KID.get_results(res, priors).ϕ_optim
    for (i, k) in enumerate(u_names)
        v = config["prior"]["parameters"][k]
        config["prior"]["parameters"][k] = (mean = u_bests[i], var = v.var, lbound = v.lbound, ubound = v.ubound)
    end
    res = KID.calibrate(KID.OptimStyle(), priors, config, ref_stats)
else
    error("Calibration method not implemented!")
end

u_bests = KID.get_results(res, priors)
KID.print_results(u_bests.ϕ_optim, u_names)
KID.save_data(
    u_bests,
    collect(u_names),
    config,
    file_name = data_save_directory * config["process"]["output_file_name"],
)
