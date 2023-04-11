import CLIMAParameters as CP
import CloudMicrophysics as CM
import Thermodynamics as TD
import Kinematic1D as KD

function get_config()
    config = Dict()
    # Define the parameter priors
    config["prior"] = get_prior_config()
    # Define parameters of observations (for validation in true-model mode)
    config["observations"] = get_observations_config()
    # Define the kalman process
    config["process"] = get_process_config()
    # Define the model
    params_calib_names = collect(keys(config["prior"]["parameters"]))
    config["model"] = get_model_config(params_calib_names)
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
    config["EKP_method"] = "UKI"
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
    config["data_names"] = ["rl", "rr"]
    # Define source of data: "file" or "perfect_model"
    config["data_source"] = "perfect_model"
    # Define number of samples for validation
    config["number_of_samples"] = 1000
    # Define random seed for generating validation samples
    config["random_seed"] = 15
    # Define the ratio of square root of covariance to G for adding artificial noise to data in the perfect-model setting
    config["scov_G_ratio"] = 0.2
    # Define offset of true values from prior means for validation
    config["true_values_offset"] = 0.25
    # Define data
    root_dir = "/Users/sajjadazimi/Postdoc/Results/01-PySDM_1D_rain_shaft/data/03-p1000/"
    config["cases"] = [(w1 = 3.0, p0 = 100000.0, Nd = 100 * 1e6, dir = root_dir * "rhow=3.0_Nd=100/")]
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

function get_model_config(params_calib_names::Array{String})
    config = Dict()
    # Define model : "KiD" or "terminal_velocity"
    config["model"] = "KiD"
    config["moisture_choice"] = "NonEquilibriumMoisture"
    config["precipitation_choice"] = "Precipitation1M"
    # Define rain formation choice: "CliMA_1M", "KK2000", "B1994", "TC1980", "LD2004", "SB2006"
    config["rain_formation_choice"] = "CliMA_1M"
    config["z_min"] = 0.0
    config["z_max"] = 3000.0
    config["n_elem"] = 64
    config["dt"] = 2.0
    config["t_ini"] = 0.0
    config["t_end"] = 3600.0
    config["dt_calib"] = 600.0
    config["t_calib"] = config["t_ini"]:config["dt_calib"]:config["t_end"]
    config["w1"] = 3.0
    config["t1"] = 600.0
    config["p0"] = 100000.0
    config["Nd"] = 100 * 1e6
    config["qtot_flux_correction"] = false
    config["r_dry"] = 0.04 * 1e-6
    config["std_dry"] = 1.4
    config["κ"] = 0.9
    config["filter"] = KD.make_filter_props(
        config["n_elem"],
        config["t_calib"];
        apply = true,
        nz_per_filtered_cell = 2,
        nt_per_filtered_cell = 120,
    )
    # Define fixed thermodynamics and microphysics parameters
    fixed_parameters = create_fixed_parameter_set(params_calib_names)
    config["thermo_params"] = fixed_parameters.thermo_params
    config["fixed_microphys_param_pairs"] = fixed_parameters.fixed_microphys_param_pairs

    return config
end

function create_fixed_parameter_set(params_calib_names::Array{String})
    FT = Float64
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    override_file = joinpath("override_dict.toml")
    open(override_file, "w") do io
        println(io, "[mean_sea_level_pressure]")
        println(io, "alias = \"MSLP\"")
        println(io, "value = 100000.0")
        println(io, "type = \"float\"")
        println(io, "[gravitational_acceleration]")
        println(io, "alias = \"grav\"")
        println(io, "value = 9.80665")
        println(io, "type = \"float\"")
        println(io, "[gas_constant]")
        println(io, "alias = \"gas_constant\"")
        println(io, "value = 8.314462618")
        println(io, "type = \"float\"")
        println(io, "[adiabatic_exponent_dry_air]")
        println(io, "alias = \"kappa_d\"")
        println(io, "value = 0.2855747338575384")
        println(io, "type = \"float\"")
        println(io, "[isobaric_specific_heat_vapor]")
        println(io, "alias = \"cp_v\"")
        println(io, "value = 1850.0")
        println(io, "type = \"float\"")
        println(io, "[molar_mass_dry_air]")
        println(io, "alias = \"molmass_dryair\"")
        println(io, "value = 0.02896998")
        println(io, "type = \"float\"")
        println(io, "[molar_mass_water]")
        println(io, "alias = \"molmass_water\"")
        println(io, "value = 0.018015")
        println(io, "type = \"float\"")
    end
    toml_dict = CP.create_toml_dict(FT; override_file, dict_type = "alias")
    isfile(override_file) && rm(override_file; force = true)

    aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TD.Parameters.ThermodynamicsParameters{FT}(; pairs...)

    aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
    aliases = setdiff(aliases, params_calib_names)
    aliases = setdiff(aliases, ["thermo_params"])
    fixed_microphys_param_pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")

    return (thermo_params = thermo_params, fixed_microphys_param_pairs = fixed_microphys_param_pairs)
end
