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
        # shared params
        # "τ_cond_evap" => (mean = 10.0, var = 1.0, lbound = 1.0, ubound = Inf),
        "a_vent_rai" => (mean = 1.5, var = 0.5, lbound = 0.0, ubound = Inf),
        # "b_vent_rai" => (mean = 0.53, var = 0.1, lbound = 0.0, ubound = Inf),
        "χv_rai" => (mean = 0.1, var = 0.03, lbound = 0.0, ubound = Inf),
        #"Δm_rai" => (mean = 0.0, var = 0.5, lbound = -3.0, ubound = 3.0),
        #"Δv_rai" => (mean = 0.0, var = 0.5, lbound = -3.0, ubound = 3.0),

        # CliMA_1M
        "q_liq_threshold" => (mean = 5.0e-4, var = 1e-4, lbound = 0.0, ubound = 2e-3),
        "k_thrshld_stpnss" => (mean = 2.0, var = 0.5, lbound = 0.5, ubound = 100.0),
        "τ_acnv_rai" => (mean = 1000.0, var = 200.0, lbound = 100.0, ubound = Inf),
        "χa_rai" => (mean = 4.0, var = 1.0, lbound = 0.0, ubound = Inf),
        #"Δa_rai" => (mean = 0.0, var = 0.5, lbound = -3.0, ubound = 3.0),

        # #TC1980
        # "k_thrshld_stpnss" => (mean = 2.0, var = 0.5, lbound = 0.5, ubound = 100.0),
        # "D_acnv_TC1980" => (mean = 3268.0, var = 100.0, lbound = 0.0, ubound = Inf), 
        # "a_acnv_TC1980" => (mean = 2.33333333333333333333, var = 0.5, lbound = 0.0, ubound = Inf), 
        # "r_0_acnv_TC1980" => (mean = 7e-6, var = 1e-6, lbound = 0.0, ubound = Inf), 
        # "A_acc_TC1980" => (mean =  4.7, var = 1.0, lbound = 0.0, ubound = Inf), 

        # # B1994
        # "k_thrshld_stpnss" => (mean = 2.0, var = 0.5, lbound = 0.5, ubound = 100.0),
        # "b_acnv_B1994" => (mean = 4.7, var = 1.0, lbound = 0.0, ubound = Inf), 
        # "d_low_acnv_B1994" => (mean = 3.9, var = 1.0, lbound = 0.0, ubound = Inf), 
        # "d_high_acnv_B1994" => (mean = 9.9, var = 1.0, lbound = 0.0, ubound = Inf), 
        # "N_0_B1994" => (mean = 2e8, var = 5e7, lbound = 0.0, ubound = Inf), 
        # "A_acc_B1994" => (mean = 6.0, var = 1.0, lbound = 0.0, ubound = Inf), 

        # # KK2000
        # "A_acnv_KK2000" => (mean = 7.42e13, var = 1e13, lbound = 0.0, ubound = Inf), 
        # "a_acnv_KK2000" => (mean = 2.47, var = 0.5, lbound = 0.0, ubound = Inf), 
        # "c_acnv_KK2000" => (mean = -1.47, var = 0.5, lbound = -Inf, ubound = Inf), 
        # "A_acc_KK2000" => (mean = 67.0, var = 10.0, lbound = 0.0, ubound = Inf), 
        # "a_acc_KK2000" => (mean = 1.15, var = 0.4, lbound = 0.0, ubound = Inf), 
        # "b_acc_KK2000" => (mean = -1.3, var = 0.3, lbound = -Inf, ubound = Inf), 

        # # LD2004
        # "k_thrshld_stpnss" => (mean = 2.0, var = 0.5, lbound = 0.5, ubound = 100.0),
        # "χa_rai" => (mean = 4.0, var = 1.0, lbound = 0.0, ubound = Inf),
        # "Δa_rai" => (mean = 0.0, var = 0.5, lbound = -3.0, ubound = 3.0),
        # "R_6C_coeff_LD2004" => (mean = 7.5, var = 1.0, lbound = -Inf, ubound = Inf),
        # "E_0_LD2004" => (mean = 1.08e10, var = 1e9, lbound = 0.0, ubound = Inf),
    )
    return config
end

function get_process_config()
    config = Dict()
    # Define method of calibration : currently only EKP and Optim are supported
    config["method"] = "EKP"
    # Define minimum number of iterations for EKP
    config["n_iter_min"] = 5
    # Define maximum number of iterations for EKP
    config["n_iter_max"] = 10
    # Define number of parameter ensemle for EKP (Inversion)
    config["n_ens"] = 15
    # Define EKP time step
    config["Δt"] = 0.5
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
    config["data_names"] = ["rlr"]
    # Define source of data: "file" or "perfect_model"
    config["data_source"] = "perfect_model"
    # Define number of samples for validation
    config["number_of_samples"] = 500
    # Define the ratio of square root of covariance to G for adding artificial noise to data in the perfect-model setting
    config["scov_G_ratio"] = 0.0
    # Define offset of true values from prior means for validation
    config["true_values_offset"] = 0.25
    # Define data
    root_dir = "/Users/sajjadazimi/Postdoc/Results/01-PySDM_1D_rain_shaft/data/03-p1000/"
    config["cases"] = [
        # (w1 = 1.0, p0 = 100000.0, Nd = 10 * 1e6, dir = root_dir * "rhow=1.0_Nd=10/"),
        # (w1 = 2.0, p0 = 100000.0, Nd = 10 * 1e6, dir = root_dir * "rhow=2.0_Nd=10/"),
        # (w1 = 3.0, p0 = 100000.0, Nd = 10 * 1e6, dir = root_dir * "rhow=3.0_Nd=10/"),
        # (w1 = 4.0, p0 = 100000.0, Nd = 10 * 1e6, dir = root_dir * "rhow=4.0_Nd=10/"),
        # (w1 = 1.0, p0 = 100000.0, Nd = 20 * 1e6, dir = root_dir * "rhow=1.0_Nd=20/"),
        # (w1 = 2.0, p0 = 100000.0, Nd = 20 * 1e6, dir = root_dir * "rhow=2.0_Nd=20/"),
        # (w1 = 3.0, p0 = 100000.0, Nd = 20 * 1e6, dir = root_dir * "rhow=3.0_Nd=20/"),
        # (w1 = 4.0, p0 = 100000.0, Nd = 20 * 1e6, dir = root_dir * "rhow=4.0_Nd=20/"),
        # (w1 = 1.0, p0 = 100000.0, Nd = 50 * 1e6, dir = root_dir * "rhow=1.0_Nd=50/"),
        # (w1 = 2.0, p0 = 100000.0, Nd = 50 * 1e6, dir = root_dir * "rhow=2.0_Nd=50/"),
        # (w1 = 3.0, p0 = 100000.0, Nd = 50 * 1e6, dir = root_dir * "rhow=3.0_Nd=50/"),
        # (w1 = 4.0, p0 = 100000.0, Nd = 50 * 1e6, dir = root_dir * "rhow=4.0_Nd=50/"),
        # (w1 = 1.0, p0 = 100000.0, Nd = 100 * 1e6, dir = root_dir * "rhow=1.0_Nd=100/"),
        # (w1 = 2.0, p0 = 100000.0, Nd = 100 * 1e6, dir = root_dir * "rhow=2.0_Nd=100/"),
        (w1 = 3.0, p0 = 100000.0, Nd = 100 * 1e6, dir = root_dir * "rhow=3.0_Nd=100/"),
        # (w1 = 4.0, p0 = 100000.0, Nd = 100 * 1e6, dir = root_dir * "rhow=4.0_Nd=100/"),
        # (w1 = 1.0, p0 = 100000.0, Nd = 200 * 1e6, dir = root_dir * "rhow=1.0_Nd=200/"),
        # (w1 = 2.0, p0 = 100000.0, Nd = 200 * 1e6, dir = root_dir * "rhow=2.0_Nd=200/"),
        # (w1 = 3.0, p0 = 100000.0, Nd = 200 * 1e6, dir = root_dir * "rhow=3.0_Nd=200/"),
        # (w1 = 4.0, p0 = 100000.0, Nd = 200 * 1e6, dir = root_dir * "rhow=4.0_Nd=200/"),
        # (w1 = 1.0, p0 = 100000.0, Nd = 500 * 1e6, dir = root_dir * "rhow=1.0_Nd=500/"),
        # (w1 = 2.0, p0 = 100000.0, Nd = 500 * 1e6, dir = root_dir * "rhow=2.0_Nd=500/"),
        # (w1 = 3.0, p0 = 100000.0, Nd = 500 * 1e6, dir = root_dir * "rhow=3.0_Nd=500/"),
        # (w1 = 4.0, p0 = 100000.0, Nd = 500 * 1e6, dir = root_dir * "rhow=4.0_Nd=500/"),
        # (w1 = 1.0, p0 = 100000.0, Nd = 1000 * 1e6, dir = root_dir * "rhow=1.0_Nd=1000/"),
        # (w1 = 2.0, p0 = 100000.0, Nd = 1000 * 1e6, dir = root_dir * "rhow=2.0_Nd=1000/"),
        # (w1 = 3.0, p0 = 100000.0, Nd = 1000 * 1e6, dir = root_dir * "rhow=3.0_Nd=1000/"),
        # (w1 = 4.0, p0 = 100000.0, Nd = 1000 * 1e6, dir = root_dir * "rhow=4.0_Nd=1000/"),
    ]
    # Define type of data
    config["data_type"] = Float64
    # Define normalization vector
    config["ynorm"] = ones(length(keys(config["cases"])), length(config["data_names"]))
    return config
end

function get_stats_config()
    config = Dict()
    # Define if pca is performed
    config["perform_pca"] = false
    # Define fraction of variance loss when performing PCA.
    config["variance_loss"] = 1e-2
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
    # Define rain formation choice: "CliMA_1M", "KK2000", "B1994", "TC1980", "LD2004"
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
        println(io, "[cloud_liquid_water_specific_humidity_autoconversion_threshold]")
        println(io, "alias = \"q_liq_threshold\"")
        println(io, "value = 5.0e-4")
        println(io, "type = \"float\"")
        println(io, "[threshold_smooth_transition_steepness]")
        println(io, "alias = \"k_thrshld_stpnss\"")
        println(io, "value = 100.0")
        println(io, "type = \"float\"")
        println(io, "[rain_autoconversion_timescale]")
        println(io, "alias = \"τ_acnv_rai\"")
        println(io, "value = 1000.0")
        println(io, "type = \"float\"")
        println(io, "[condensation_evaporation_timescale]")
        println(io, "alias = \"τ_cond_evap\"")
        println(io, "value = 10.0")
        println(io, "type = \"float\"")
        println(io, "[rain_ventillation_coefficient_a]")
        println(io, "alias = \"a_vent_rai\"")
        println(io, "value = 1.5")
        println(io, "type = \"float\"")
        println(io, "[rain_ventillation_coefficient_b]")
        println(io, "alias = \"b_vent_rai\"")
        println(io, "value = 0.53")
        println(io, "type = \"float\"")
        println(io, "[rain_mass_size_relation_coefficient_chim]")
        println(io, "alias = \"χm_rai\"")
        println(io, "value = 1.0")
        println(io, "type = \"float\"")
        println(io, "[rain_cross_section_size_relation_coefficient_chia]")
        println(io, "alias = \"χa_rai\"")
        println(io, "value = 1.0")
        println(io, "type = \"float\"")
        println(io, "[rain_terminal_velocity_size_relation_coefficient_chiv]")
        println(io, "alias = \"χv_rai\"")
        println(io, "value = 1.0")
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
