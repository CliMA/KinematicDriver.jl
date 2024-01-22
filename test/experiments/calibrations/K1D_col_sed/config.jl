import CLIMAParameters as CP
import CloudMicrophysics as CM
import Thermodynamics as TD
import Kinematic1D.CalibrateCMP as KCP

function get_config()
    config = Dict()
    # Define the parameter priors
    config["prior"] = get_prior_config()
    # Define parameters of observations (for validation in true-model mode)
    config["observations"] = get_observations_config()
    # Define the kalman process
    config["process"] = get_process_config()
    # Define the model
    config["model"] = get_model_config()
    # Define statistics
    config["statistics"] = get_stats_config()
    return config
end

function get_prior_config()
    config = Dict()
    # Define prior mean and bounds on the parameters.
    config["parameters"] = Dict(
        "χv_rai" => (mean = 0.1, var = 0.03, lbound = 0.0, ubound = 1.0),
        "χa_rai" => (mean = 4.0, var = 1.0, lbound = 0.0, ubound = 10.0),
    )
    return config
end

function get_process_config()
    config = Dict()
    # Define method of calibration : currently only EKP and Optim are supported
    config["method"] = "EKP"
    # Define mini batch size for EKP
    config["batch_size"] = 1
    # Define number of iterations for EKP
    config["n_iter"] = 10
    # Define number of parameter ensemle for EKP (Inversion)
    config["n_ens"] = 15
    # Define EKP time step
    config["Δt"] = 1.0
    config["EKP_method"] = "EKI"
    # Choose regularization factor α ∈ (0,1] for UKI, when enough observation data α=1: no regularization
    config["α_reg"] = 1.0
    # UKI parameter
    # update_freq = 1 : approximate posterior covariance matrix with an uninformative prior
    #               0 : weighted average between posterior covariance matrix with an uninformative prior and prior
    config["update_freq"] = 1
    # Define Optim absolute tolerance for convergence
    config["tol"] = 1e-3
    # Define Optim maximum iterations
    config["maxiter"] = 20000
    # Define output file name
    config["output_file_name"] = "parameters_EKP.jld2"
    return config
end

function get_observations_config()
    config = Dict()
    # Define data names.
    config["data_names"] = ["rl", "rr", "Nl", "Nr"]
    # Define source of data: "file" or "perfect_model"
    config["data_source"] = "file"
    # Define number of samples for validation
    config["number_of_samples"] = 1000
    # Define random seed for generating validation samples
    config["random_seed"] = 15
    # Define the ratio of square root of covariance to G for adding artificial noise to data in the perfect-model setting
    config["scov_G_ratio"] = 0.2
    # Define offset of true values from prior means for validation
    config["true_values_offset"] = 0.25
    # Define data
    root_dir = "/Users/sajjadazimi/Postdoc/Results/05-PySDM_1D_rain_shaft_col_sed/data/01-Geometric/"
    config["cases"] = [(qt = 1.6e-3, Nd = 80 * 1e6, k = 2.0, dir = root_dir * "qt=1.6_Nd=80_k=2/")]
    # Define type of data
    config["data_type"] = Float64
    return config
end

function get_stats_config()
    config = Dict()
    # Define normalization method: mean_normalized or std_normalized
    config["normalization"] = "std_normalized"
    # Define if pca is performed
    config["perform_pca"] = true
    # Define fraction of variance loss when performing PCA.
    config["variance_loss"] = 0.01
    # Define tikhonov mode: absulute or relative
    config["tikhonov_mode"] = "absolute"
    # Define tikhonov noise
    config["tikhonov_noise"] = config["perform_pca"] ? 0.0 : 1e-3
    return config
end

function get_model_config()
    config = Dict()
    config["model"] = "KiD_col_sed"
    config["qt"] = 1.6e-3
    config["Nd"] = 80 * 1e6
    config["k"] = 2.0
    config["rhod"] = 1.0
    config["precipitation_choice"] = "Precipitation2M"
    # Define rain formation choice: "CliMA_1M", "KK2000", "B1994", "TC1980", "LD2004", "VarTimeScaleAcnv", "SB2006"
    config["rain_formation_choice"] = "SB2006"
    # Define sedimentation choice: "CliMA_1M", "Chen2022", "SB2006"
    config["sedimentation_choice"] = "SB2006"
    config["precip_sources"] = true
    config["precip_sinks"] = false
    config["z_min"] = 0.0
    config["z_max"] = 3000.0
    config["n_elem"] = 64
    config["dt"] = 2.0
    config["t_ini"] = 0.0
    config["t_end"] = 3600.0
    config["dt_calib"] = 150.0
    config["t_calib"] = config["t_ini"]:config["dt_calib"]:config["t_end"]
    config["filter"] = KCP.make_filter_props(
        config["n_elem"],
        config["t_calib"];
        apply = true,
        nz_per_filtered_cell = 4,
        nt_per_filtered_cell = 30,
    )
    # Define default parameters
    FT = Float64
    config["toml_dict"] = CP.create_toml_dict(FT, dict_type = "alias")

    return config
end
