using Plots

import KinematicDriver.CalibrateCMP as KCP

include("./config.jl")

data_save_directory = KCP.make_output_directories(joinpath(@__DIR__, "run_model_output"))
output_file_name_base = "box_model_vs_obs"

config = get_config()
u_names = collect(keys(config["prior"]["parameters"]))
u_values = collect([v.mean for v in values(config["prior"]["parameters"])])

obs = KCP.get_obs!(config)
config["statistics"]["normalization"] = "mean_normalized"
ref_stats_list = KCP.make_ref_stats_list(obs, config["statistics"], KCP.get_numbers_from_config(config)...)
ref_stats = KCP.combine_ref_stats(ref_stats_list)
G = KCP.run_dyn_model(u_values, u_names, config)
Gn = KCP.normalize_sim(G, ref_stats)

model_error = KCP.compute_error_metrics(u_values, u_names, config, ref_stats)
println("loss = ", model_error.loss, ",\t mse_m = ", model_error.mse_m, ",\t mse_s = ", model_error.mse_s)

KCP.plot_box_results(
    Gn,
    ref_stats.y_full,
    sqrt.(diag(ref_stats.Γ_full)),
    config,
    path = data_save_directory,
    file_base = output_file_name_base,
)
plot!()
