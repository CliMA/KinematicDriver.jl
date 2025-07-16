import Plots
include("../Calibration.jl")

config_dir = "K1D"
data_save_directory = make_output_directories(joinpath(@__DIR__, config_dir, "run_model_output"))
output_file_name_base = "model_vs_obs_contours"

include(joinpath(@__DIR__, config_dir, "config.jl"))
config = get_config()
u_names = collect(keys(config["prior"]["parameters"]))
u_values = collect([v.mean for v in values(config["prior"]["parameters"])])

obs = get_obs!(config)
config["statistics"]["normalization"] = "mean_normalized"
ref_stats_list = make_ref_stats_list(obs, config["statistics"], get_numbers_from_config(config)...)
ref_stats = combine_ref_stats(ref_stats_list)
G = run_dyn_model(u_values, u_names, config)
Gn = normalize_sim(G, ref_stats)

model_error = compute_error_metrics(u_values, u_names, config, ref_stats)
println("loss = ", model_error.loss, ",\t mse_m = ", model_error.mse_m, ",\t mse_s = ", model_error.mse_s)

@info "Plotting"
if config["model"]["model"] == "Box"
    plot_box_results(
        Gn,
        ref_stats.y_full,
        sqrt.(diag(ref_stats.Γ_full)),
        config,
        path = data_save_directory,
        file_base = output_file_name_base,
    )
else
    compare_model_and_obs_contours(
        Gn,
        ref_stats.y_full,
        config,
        levels = range(minimum([0, minimum(Gn)]), maximum([1, maximum(Gn)]), 10),
        linewidth = 0.5,
        G_title = "KiD",
        obs_title = "PySDM",
        path = data_save_directory,
        file_base = output_file_name_base,
    )
end

Plots.plot(ref_stats.y_full)
Plots.plot!(Gn, legend = false)
