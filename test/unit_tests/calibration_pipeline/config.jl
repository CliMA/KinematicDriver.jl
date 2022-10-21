import CLIMAParameters as CP
import CloudMicrophysics as CM
import Thermodynamics as TD

function get_config()
    config = Dict()
    config["prior"] = get_prior_config()
    config["observations"] = get_observations_config()
    config["process"] = get_process_config()
    params_calib_names = collect(keys(config["prior"]["parameters"]))
    config["model"] = get_model_config(params_calib_names)
    config["statistics"] = get_stats_config()
    return config
end

function get_prior_config()
    config = Dict()
    config["parameters"] = Dict(
        "q_liq_threshold" => (mean = 5.0e-4, var = 1e-4, lbound = 0.0, ubound = 2e-3),
        "τ_acnv_rai" => (mean = 1000.0, var = 200.0, lbound = 100.0, ubound = Inf),
    )
    return config
end

function get_process_config()
    config = Dict()
    config["n_iter_min"] = 10
    config["n_iter_max"] = 20
    config["n_ens"] = 15
    config["Δt"] = 0.5
    config["tol"] = 1e-8
    config["maxiter"] = 1000
    return config
end

function get_observations_config()
    config = Dict()
    config["data_names"] = ["rlr", "qt"]
    config["data_source"] = "perfect_model"
    config["number_of_samples"] = 50
    config["scov_G_ratio"] = 0.0
    config["true_values_offset"] = 0.25
    config["cases"] = [
        (w1 = 2.0, p0 = 99000.0, Nd = 50 * 1e6, dir = "./case1/"),
        (w1 = 3.0, p0 = 99000.0, Nd = 50 * 1e6, dir = "./case2/"),
    ]
    config["data_type"] = Float64
    config["ynorm"] = ones(length(keys(config["cases"])), length(config["data_names"]))
    return config
end

function get_stats_config()
    config = Dict()
    config["perform_pca"] = false
    config["variance_loss"] = 1e-2
    config["tikhonov_mode"] = "absolute"
    config["tikhonov_noise"] = config["perform_pca"] ? 0.0 : 1e-3
    return config
end

function get_model_config(params_calib_names::Array{String})
    config = Dict()
    config["model"] = "KiD"
    config["moisture_choice"] = "NonEquilibriumMoisture"
    config["precipitation_choice"] = "Precipitation1M"
    config["rain_formation_choice"] = "CliMA_1M"
    config["z_min"] = 0.0
    config["z_max"] = 3000.0
    config["n_elem"] = 16
    config["dt"] = 4.0
    config["t_ini"] = 0.0
    config["t_end"] = 60.0
    config["dt_output"] = 60.0
    config["t_calib"] = 0:config["dt_output"]:config["t_end"]
    config["w1"] = 3.0
    config["t1"] = 600.0
    config["p0"] = 99000.0
    config["Nd"] = 50 * 1e6
    config["qtot_flux_correction"] = false
    fixed_parameters = create_fixed_parameter_set(params_calib_names)
    config["thermo_params"] = fixed_parameters.thermo_params
    config["fixed_microphys_param_pairs"] = fixed_parameters.fixed_microphys_param_pairs

    return config
end

function create_fixed_parameter_set(params_calib_names::Array{String})
    FT = Float64
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TD.Parameters.ThermodynamicsParameters{FT}(; pairs...)

    aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
    aliases = setdiff(aliases, params_calib_names)
    aliases = setdiff(aliases, ["thermo_params"])
    fixed_microphys_param_pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")

    return (thermo_params = thermo_params, fixed_microphys_param_pairs = fixed_microphys_param_pairs)
end
