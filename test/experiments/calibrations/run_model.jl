using Plots

import Kinematic1D
const KID = Kinematic1D

include("./config.jl")

data_save_directory = KID.make_output_directories()
output_file_name_base = "model_vs_obs_contours"

config = get_config()
u_names = collect(keys(config["prior"]["parameters"]))
u_values = collect([v.mean for v in values(config["prior"]["parameters"])])

truth = KID.get_obs!(config)
G = KID.run_dyn_model(u_values, u_names, config)

KID.compare_model_and_obs_contours(
    G[:],
    truth.mean,
    config,
    levels = range(minimum([0, minimum(G)]), maximum([1, maximum(G)]), 10),
    linewidth = 0.5,
    overlay = false,
    G_title = "KiD",
    obs_title = "PySDM",
    path = data_save_directory,
    file_base = output_file_name_base,
)

ref_stats = KID.ReferenceStatistics(truth, config["statistics"])
println("loss: ", KID.compute_loss(u_values, u_names, config, ref_stats))

plot(truth.mean)
plot!(G, legend = false)
