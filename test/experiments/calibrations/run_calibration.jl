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
else
    error("Calibration method not implemented!")
end

u_bests = KID.get_results(res, priors)
KID.print_results(u_bests, u_names)
KID.save_data(u_bests, config, file_name = data_save_directory * config["process"]["output_file_name"])
