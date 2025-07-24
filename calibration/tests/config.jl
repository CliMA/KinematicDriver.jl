import ClimaParams as CP
import CloudMicrophysics as CM
import Thermodynamics as TD
include("../HelperFuncs.jl")

function get_config()
    config = Dict()
    config["prior"] = get_prior_config()
    config["observations"] = get_observations_config()
    config["process"] = get_process_config()
    config["model"] = get_model_config()
    config["statistics"] = get_stats_config()
    return config
end

function get_prior_config()
    config = Dict()
    config["parameters"] = Dict(
        "rain_terminal_velocity_size_relation_coefficient_chiv" =>
            (mean = 1.0, var = 0.25, lbound = 0.0, ubound = 2.0),
        "rain_cross_section_size_relation_coefficient_chia" => (mean = 1.0, var = 0.25, lbound = 0.0, ubound = Inf),
    )
    return config
end

function get_process_config()
    config = Dict()
    config["batch_size"] = 1
    config["n_iter"] = 5
    config["n_ens"] = 10
    config["Δt"] = 1.0
    config["EKP_method"] = "EKI"
    config["augmented"] = false
    config["α_reg"] = 1.0
    config["update_freq"] = 1
    config["tol"] = 1e-8
    config["maxiter"] = 1000
    return config
end

function get_observations_config()
    config = Dict()
    config["data_names"] = ["rl", "qr"]
    config["data_source"] = "perfect_model"
    config["number_of_samples"] = 100
    config["random_seed"] = 15
    config["scov_G_ratio"] = 0.2
    config["true_values_offset"] = 0.25
    config["cases"] = [
        (w1 = 2.0, p0 = 99000.0, Nd = 50 * 1e6, dir = "./case1/"),
        (w1 = 3.0, p0 = 99000.0, Nd = 50 * 1e6, dir = "./case2/"),
    ]
    config["data_type"] = Float64
    return config
end

function get_stats_config()
    config = Dict()
    config["normalization"] = "std_normalized"
    config["perform_pca"] = false
    config["variance_loss"] = 1e-2
    config["tikhonov_mode"] = "absolute"
    config["tikhonov_noise"] = config["perform_pca"] ? 0.0 : 1e-3
    return config
end

function get_model_config()
    config = Dict()
    config["model"] = "KiD"
    config["moisture_choice"] = "NonEquilibriumMoisture"
    config["precipitation_choice"] = "Precipitation1M"
    config["rain_formation_choice"] = "CliMA_1M"
    config["sedimentation_choice"] = "CliMA_1M"
    config["precip_sources"] = true
    config["precip_sinks"] = true
    config["z_min"] = 0.0
    config["z_max"] = 3000.0
    config["n_elem"] = 10
    config["dt"] = 10.0
    config["t_ini"] = 0.0
    config["t_end"] = 1200.0
    config["dt_calib"] = 300.0
    config["t_calib"] = config["t_ini"]:config["dt_calib"]:config["t_end"]
    config["w1"] = 3.0
    config["t1"] = 600.0
    config["p0"] = 99000.0
    config["Nd"] = 50 * 1e6
    config["qtot_flux_correction"] = false
    config["open_system_activation"] = false
    config["local_activation"] = false
    config["r_dry"] = 0.04 * 1e-6
    config["std_dry"] = 1.4
    config["κ"] = 0.9
    config["filter"] = make_filter_props([config["n_elem"], config["n_elem"]], config["t_calib"]; apply = false)
    # Define default parameters
    params = create_parameter_set()
    config["toml_dict"] = params.toml_dict
    config["thermo_params"] = params.thermo_params
    config["air_params"] = params.air_params
    config["activation_params"] = params.activation_params

    return config
end

function create_parameter_set()
    FT = Float64
    toml_dict = CP.create_toml_dict(FT)
    thermo_params = TD.Parameters.ThermodynamicsParameters(toml_dict)
    air_params = CM.Parameters.AirProperties(toml_dict)
    activation_params = CM.Parameters.AerosolActivationParameters(toml_dict)
    return (; toml_dict, thermo_params, air_params, activation_params)
end
