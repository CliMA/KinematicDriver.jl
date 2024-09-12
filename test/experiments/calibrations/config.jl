import ClimaParams as CP
import CloudMicrophysics as CM
import Thermodynamics as TD
import KinematicDriver.CalibrateCMP as KCP

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
        "SB2006_collection_kernel_coeff_kcc" =>
            (mean = 4.44 * 1e9, var = 0.3 * 8.88 * 1e8, lbound = 0.0, ubound = Inf),
        "SB2006_collection_kernel_coeff_kcr" =>
           (mean = 5.25, var = 0.35 * 1.05, lbound = 0.0, ubound = Inf),
        "SB2006_collection_kernel_coeff_krr" => 
            (mean = 7.12, var = 0.3 * 1.424, lbound = 0.0, ubound = Inf),
        "SB2006_raindrops_terminal_velocity_coeff_aR" => 
           (mean = 9.65, var = 0.25, lbound = 0.0, ubound = 10.2),
    )
    return config
end

function get_process_config()
    config = Dict()
    # Define method of calibration : currently only EKP and Optim are supported
    config["method"] = "EKP"
    # Define mini batch size for EKP
    config["batch_size"] = 13
    # Define number of iterations for EKP
    config["n_iter"] = 15
    # Define number of parameter ensemble for EKP (Inversion)
    config["n_ens"] = 20
    # Define EKP time step
    config["Δt"] = 0.1
    config["EKP_method"] = "ETKI"
    # Define whether state vector is augmented by parameters for Bayesian regularization
    config["augmented"] = false
    # Define Bayesian regularization degree scale (prior cov * scale)
    config["reg_degree_scale"] = 1.0
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
    config["data_names"] = ["reff_top", "Z", "rainrate_surface"]
    # Define source of data: "file" or "perfect_model"
    config["data_source"] = "file"
    # Define number of samples for validation
    config["number_of_samples"] = 1000
    # Define random seed for generating validation samples
    config["random_seed"] = 15
    # Define the ratio of square root of covariance to G for adding artificial noise to data in the perfect-model setting
    config["scov_G_ratio"] = 0.2
    # Define offset of true values from prior means for validation
    config["true_values_offset"] = 0.1
    # Define data
    root_dir = "/Users/caterinacroci/Desktop/data/"
    config["cases"] =[  # w = 2.0 m/s
                        #=(w1 = 2.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.1,
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=990_Nd=50/"),
                        (w1 = 2.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.1,
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=990_Nd=100/"),
                        (w1 = 2.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=990_Nd=500/"),
                        (w1 = 2.0, p0 = 100000.0, Nd = 50 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1000_Nd=50/"),
                        (w1 = 2.0, p0 = 100000.0, Nd = 100 * 1e6, std_dry = 1.1,  
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1000_Nd=100/"),
                        (w1 = 2.0, p0 = 100000.0, Nd = 500 * 1e6, std_dry = 1.1,  
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1000_Nd=500/"),
                        (w1 = 2.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1007_Nd=50/"),
                        (w1 = 2.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1007_Nd=100/"),

                        (w1 = 2.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=990_Nd=50/"),
                        (w1 = 2.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=990_Nd=100/"),
                        (w1 = 2.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=990_Nd=500/"),
                        (w1 = 2.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.8,  
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=1007_Nd=50/"),
                        (w1 = 2.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=1007_Nd=100/"),
                        
                        # w = 4.0 m/s
                        (w1 = 4.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=990_Nd=50/"),
                        (w1 = 4.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=990_Nd=100/"),
                        (w1 = 4.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=990_Nd=500/"),
                        (w1 = 4.0, p0 = 100000.0, Nd = 50 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1000_Nd=50/"),
                        (w1 = 4.0, p0 = 100000.0, Nd = 100 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1000_Nd=100/"),
                        (w1 = 4.0, p0 = 100000.0, Nd = 500 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1000_Nd=500/"),
                        (w1 = 4.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1007_Nd=50/"),
                        (w1 = 4.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1007_Nd=100/"),

                        (w1 = 4.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=990_Nd=50/"),
                        (w1 = 4.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=990_Nd=100/"),
                        (w1 = 4.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=990_Nd=500/"),
                        (w1 = 4.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=1007_Nd=50/"),
                        (w1 = 4.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=1007_Nd=100/"),

                        # w = 5.0 m/s
                        (w1 = 5.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=990_Nd=50/"),
                        (w1 = 5.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.1,
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=990_Nd=100/"),
                        (w1 = 5.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.1,
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=990_Nd=500/"),
                        (w1 = 5.0, p0 = 100000.0, Nd = 50 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1000_Nd=50/"),
                        (w1 = 5.0, p0 = 100000.0, Nd = 100 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1000_Nd=100/"),
                        (w1 = 5.0, p0 = 100000.0, Nd = 500 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1000_Nd=500/"),
                        (w1 = 5.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1007_Nd=50/"),
                        (w1 = 5.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.1, 
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1007_Nd=100/"),

                        (w1 = 5.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=990_Nd=50/"),
                        (w1 = 5.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=990_Nd=100/"),
                        (w1 = 5.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=990_Nd=500/"),
                        (w1 = 5.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.8,
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=1007_Nd=50/"),
                        (w1 = 5.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.8, 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=1007_Nd=100/"),=#

                        (w1 = 2.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.1, t_cal = [620.0, 1065.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=990_Nd=50/"),
                        (w1 = 2.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.1, t_cal = [835.0, 1250.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=990_Nd=100/"),
                        (w1 = 2.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.1, t_cal = [1310.0, 1760.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=990_Nd=500/"),
                        (w1 = 2.0, p0 = 100000.0, Nd = 50 * 1e6, std_dry = 1.1, t_cal = [920.0, 1697.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1000_Nd=50/"),
                        (w1 = 2.0, p0 = 100000.0, Nd = 100 * 1e6, std_dry = 1.1,  t_cal = [1295.0, 1900.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1000_Nd=100/"),
                        (w1 = 2.0, p0 = 100000.0, Nd = 500 * 1e6, std_dry = 1.1,  t_cal = [1960.0, 2655.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1000_Nd=500/"),
                        (w1 = 2.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.1, t_cal = [1250.0, 2350.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1007_Nd=50/"),
                        (w1 = 2.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.1, t_cal = [1500.0, 2700.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_2/mean=0.04_std=1.1_p0=1007_Nd=100/"),

                        (w1 = 2.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.8, t_cal = [600.0, 1150.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=990_Nd=50/"),
                        (w1 = 2.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.8, t_cal = [835.0, 1375.0], 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=990_Nd=100/"),
                        (w1 = 2.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.8, t_cal = [1310.0, 1775.0], 
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=990_Nd=500/"),
                        (w1 = 2.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.8, t_cal = [1000.0, 2070.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=1007_Nd=50/"),
                        (w1 = 2.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.8, t_cal = [1500.0, 2575.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_2/mean=0.06_std=1.8_p0=1007_Nd=100/"),

                        # w = 4.0 m/s
                        (w1 = 4.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.1, t_cal = [620.0, 1025.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=990_Nd=50/"),
                        (w1 = 4.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.1, t_cal = [620.0, 1125.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=990_Nd=100/"),
                        (w1 = 4.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.1, t_cal = [1020.0, 1480.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=990_Nd=500/"),
                        (w1 = 4.0, p0 = 100000.0, Nd = 50 * 1e6, std_dry = 1.1, t_cal = [600.0, 1400.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1000_Nd=50/"),
                        (w1 = 4.0, p0 = 100000.0, Nd = 100 * 1e6, std_dry = 1.1, t_cal = [880.0, 1415.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1000_Nd=100/"),
                        (w1 = 4.0, p0 = 100000.0, Nd = 500 * 1e6, std_dry = 1.1, t_cal = [1440.0, 2030.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1000_Nd=500/"),
                        (w1 = 4.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.1, t_cal = [825.0, 1782.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1007_Nd=50/"),
                        (w1 = 4.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.1, t_cal = [1205.0, 2030.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_4/mean=0.04_std=1.1_p0=1007_Nd=100/"),

                        (w1 = 4.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.8, t_cal = [620.0, 1100.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=990_Nd=50/"),
                        (w1 = 4.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.8, t_cal = [620.0, 1120.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=990_Nd=100/"),
                        (w1 = 4.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.8, t_cal = [1000.0, 1460.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=990_Nd=500/"),
                        (w1 = 4.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.8, t_cal = [825.0, 1800.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=1007_Nd=50/"),
                        (w1 = 4.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.8, t_cal = [1205.0, 1925.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_4/mean=0.06_std=1.8_p0=1007_Nd=100/"),

                        # w = 5.0 m/s
                        (w1 = 5.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.1, t_cal = [620.0, 1050.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=990_Nd=50/"),
                        (w1 = 5.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.1, t_cal = [620.0, 1200.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=990_Nd=100/"),
                        (w1 = 5.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.1, t_cal = [925.0, 1407.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=990_Nd=500/"),
                        (w1 = 5.0, p0 = 100000.0, Nd = 50 * 1e6, std_dry = 1.1, t_cal = [620.0, 1200.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1000_Nd=50/"),
                        (w1 = 5.0, p0 = 100000.0, Nd = 100 * 1e6, std_dry = 1.1, t_cal = [755.0, 1415.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1000_Nd=100/"),
                        (w1 = 5.0, p0 = 100000.0, Nd = 500 * 1e6, std_dry = 1.1, t_cal = [1240.0, 1825.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1000_Nd=500/"),
                        (w1 = 5.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.1, t_cal = [620.0, 1675.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1007_Nd=50/"),
                        (w1 = 5.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.1, t_cal = [1020.0, 1812.0],
                            r_dry = 0.04 * 1e-6, dir = root_dir * "rhow_5/mean=0.04_std=1.1_p0=1007_Nd=100/"),

                        (w1 = 5.0, p0 = 99000.0, Nd = 50 * 1e6, std_dry = 1.8, t_cal = [620.0, 1000.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=990_Nd=50/"),
                        (w1 = 5.0, p0 = 99000.0, Nd = 100 * 1e6, std_dry = 1.8, t_cal = [620.0, 1025.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=990_Nd=100/"),
                        (w1 = 5.0, p0 = 99000.0, Nd = 500 * 1e6, std_dry = 1.8, t_cal = [850.0, 1300.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=990_Nd=500/"),
                        (w1 = 5.0, p0 = 100700.0, Nd = 50 * 1e6, std_dry = 1.8, t_cal = [620.0, 1550.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=1007_Nd=50/"),
                        (w1 = 5.0, p0 = 100700.0, Nd = 100 * 1e6, std_dry = 1.8, t_cal = [1020.0, 1680.0],
                            r_dry = 0.06 * 1e-6, dir = root_dir * "rhow_5/mean=0.06_std=1.8_p0=1007_Nd=100/"),]
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
    config["model"] = "KiD"
    config["moisture_choice"] = "NonEquilibriumMoisture"
    config["precipitation_choice"] = "Precipitation2M"
    # Define rain formation choice: "CliMA_1M", "KK2000", "B1994", "TC1980", "LD2004", "VarTimeScaleAcnv", "SB2006", "SB2006NL"
    config["rain_formation_choice"] = "SB2006NL"
    # Define sedimentation choice: "CliMA_1M", "Chen2022", "SB2006"
    config["sedimentation_choice"] = "SB2006"
    config["precip_sources"] = true 
    config["precip_sinks"] = true
    config["z_min"] = 0.0
    config["z_max"] = 4000.0
    config["n_elem"] = 80
    config["dt"] = 0.75
    config["t_ini"] = 0.0
    config["t_end"] = 3600.0
    config["dt_calib"] = 120.0
    config["t_calib"] = 600.0:config["dt_calib"]:config["t_end"]
    config["w1"] = 3.0
    config["t1"] = 600.0
    config["p0"] = 99000.0
    config["Nd"] = 50 * 1e6
    config["qtot_flux_correction"] = false
    config["open_system_activation"] = true
    config["r_dry"] = 0.04 * 1e-6
    config["std_dry"] = 1.1
    config["κ"] = 0.9
    config["filter"] = KCP.make_filter_props(
        [1, config["n_elem"], 1], # nz (for each variable)
        config["t_calib"];
        apply = true,
        nz_per_filtered_cell = [1, 10, 1],
        nt_per_filtered_cell = 120,
    )
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
        println(io, "[SB2006_raindrops_min_mass]")
        println(io, "alias = \"raindrops_min_mass\"")
        println(io, "value = 6.54e-11")
        println(io, "[alpha]")
        println(io, "alias = \"alpha\"")
        println(io, "value = 1.0")
        println(io, "type = \"float\"")
    end
    toml_dict = CP.create_toml_dict(FT; override_file)
    isfile(override_file) && rm(override_file; force = true)

    thermo_params = TD.Parameters.ThermodynamicsParameters(toml_dict)
    air_params = CM.Parameters.AirProperties(toml_dict)
    activation_params = CM.Parameters.AerosolActivationParameters(toml_dict)

    return (; toml_dict, thermo_params, air_params, activation_params)
end
