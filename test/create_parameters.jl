import KinematicDriver.Common as CO
import KinematicDriver.K1DModel as KID
import ClimaParams as CP
import CloudMicrophysics as CM
import Thermodynamics as TD
import Cloudy as CL

#! format: off
function override_toml_dict(
    out_dir::String,
    toml_dict::CP.AbstractTOMLDict;
    w1 = 2.0,
    t1 = 600.0,
    p0 = 100700.0,
    z_0 = 0.0,
    z_1 = 740.0,
    z_2 = 3260.0,
    rv_0 = 0.015,
    rv_1 = 0.0138,
    rv_2 = 0.0024,
    tht_0 = 297.9,
    tht_1 = 297.9,
    tht_2 = 312.66,
    precip_sources = 1,
    precip_sinks = 1,
    qtot_flux_correction = 0,
    prescribed_Nd = 100 * 1e6,
    r_dry = 0.04 * 1e-6,
    std_dry = 1.4,
    κ = 1.12,
)
    FT = CP.float_type(toml_dict)
    override_file = Dict(
        "mean_sea_level_pressure" => Dict("value" => 100000.0, "type" => "float"),
        "gravitational_acceleration" => Dict("value" => 9.80665, "type" => "float"),
        "gas_constant" => Dict("value" => 8.314462618, "type" => "float"),
        "adiabatic_exponent_dry_air" => Dict("value" => 0.2855747338575384, "type" => "float"),
        "isobaric_specific_heat_vapor" => Dict("value" => 1850.0, "type" => "float"),
        "molar_mass_dry_air" => Dict("value" => 0.02896998, "type" => "float"),
        "molar_mass_water" => Dict("value" => 0.018015, "type" => "float"),
        "cloud_liquid_water_specific_humidity_autoconversion_threshold" => Dict("value" => 0.0001, "type" => "float"),
        "prescribed_flow_w1" => Dict("value" => 1.0, "type" => "float"),
        "prescribed_flow_t1" => Dict("value" => t1, "type" => "float"),
        "surface_pressure" => Dict("value" => 99000.0, "type" => "float"),
        "precipitation_sources_flag" => Dict("value" => precip_sources, "type" => "bool"),
        "precipitation_sinks_flag" => Dict("value" => precip_sinks, "type" => "bool"),
        "qtot_flux_correction_flag" => Dict("value" => qtot_flux_correction, "type" => "bool"),
        "prescribed_Nd" => Dict("value" => 50 * 1e6, "type" => "float"),
        "r_dry" => Dict("value" => r_dry, "type" => "float"),
        "std_dry" => Dict("value" => 1.1, "type" => "float"),
        "kappa" => Dict("value" => 0.9, "type" => "float"),
        "init_cond_z0" => Dict("value" => z_0, "type" => "float"),
        "init_cond_z1" => Dict("value" => z_1, "type" => "float"),
        "init_cond_z2" => Dict("value" => z_2, "type" => "float"),
        "init_cond_rv0" => Dict("value" => rv_0, "type" => "float"),
        "init_cond_rv1" => Dict("value" => rv_1, "type" => "float"),
        "init_cond_rv2" => Dict("value" => rv_2, "type" => "float"),
        "init_cond_theta0" => Dict("value" => tht_0, "type" => "float"),
        "init_cond_theta1" => Dict("value" => tht_1, "type" => "float"),
        "init_cond_theta2" => Dict("value" => tht_2, "type" => "float"),
        "SB2006_raindrops_min_mass" => Dict("value" => 6.54e-11, "type" => "float"),
        "SB2006_raindrops_size_distribution_coeff_N0_min" => Dict("value" => 3.5e5, "type" => "float"),
        "SB2006_raindrops_size_distribution_coeff_N0_max" => Dict("value" => 2e11, "type" => "float"),
        "SB2006_raindrops_size_distribution_coeff_lambda_min" => Dict("value" => 1e3, "type" => "float"),
        "SB2006_raindrops_size_distribution_coeff_lambda_max" => Dict("value" => 4e4, "type" => "float"),
        # NEW!!!!
        #"SB2006_cloud_gamma_distribution_coeff_mu" => Dict("value" => 0.3578, "type" => "float"),
        #"SB2006_collection_kernel_coeff_krr" => Dict("value" => 9.657, "type" => "float"),
        #"SB2006_collection_kernel_coeff_kcc" => Dict("value" => 6.719e9, "type" => "float"),
        #"SB2006_collection_kernel_coeff_kcr" => Dict("value" => 3.998, "type" => "float"),
        #"SB2006_autoconversion_correcting_function_coeff_A" => Dict("value" => 360.1, "type" => "float"),
        #"SB2006_cloud_gamma_distribution_parameter" => Dict("value" => 5.0, "type" => "float"),
        #"SB2006_autoconversion_correcting_function_coeff_a" => Dict("value" => 0.3395, "type" => "float"),
        #"SB2006_autoconversion_correcting_function_coeff_b" => Dict("value" => 8.786, "type" => "float"),
        #"SB2006_accretion_correcting_function_coeff_c" => Dict("value" => 2, "type" => "float"),
        #"SB2006_collection_kernel_coeff_kapparr" => Dict("value" => 5.566, "type" => "float"),
        #"SB2006_raindrops_terminal_velocity_coeff_aR" => Dict("value" => 9.208, "type" => "float"),
        #"SB2006_raindrops_terminal_velocity_coeff_bR" => Dict("value" => 10.32, "type" => "float"),
        #"SB2006_raindrops_terminal_velocity_coeff_cR" => Dict("value" => 304.8, "type" => "float"),    
    )
    toml_dict = CP.create_toml_dict(FT; override_file)
    return toml_dict
end

function create_thermodynamics_parameters(toml_dict)
    return TD.Parameters.ThermodynamicsParameters(toml_dict)
end

function create_common_parameters(toml_dict)
    common_params = CO.Parameters.CommonParameters(toml_dict)
    if !isbits(common_params)
        print(common_params)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return common_params
end

function create_kid_parameters(toml_dict)
    kid_params = KID.Parameters.KinematicDriverParameters(toml_dict)
    if !isbits(kid_params)
        print(kid_params)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return kid_params
end

function create_cloudy_parameters(FT, dist_names::NTuple{ND, String} = ("gamma", "gamma")) where {ND}

    # Create water category distributions
    function make_dist(dist_name)
        if dist_name == "monodisperse"
            return CL.ParticleDistributions.MonodispersePrimitiveParticleDistribution(FT(0), FT(1))
        elseif dist_name == "exponential"
            return CL.ParticleDistributions.GammaPrimitiveParticleDistribution(FT(0), FT(1))
        elseif dist_name == "gamma"
            return CL.ParticleDistributions.GammaPrimitiveParticleDistribution(FT(0), FT(1), FT(1))
        end
    end
    pdists::NTuple{ND, CL.ParticleDistributions.PrimitiveParticleDistribution} = map(dist_names) do name
        make_dist(name)
    end

    # Create tuple of number of prognostic moments for each distribution
    NProgMoms::NTuple{ND, Int} = map(pdists) do dist
        CL.ParticleDistributions.nparams(dist)
    end

    # Define norms of number per volume and mass of particles
    norms::Tuple{FT, FT} = (FT(1e6), FT(1e-9)) # 1e6 1/m^3, 1e-9 kg

    # Compute normalizing factors of prognostic moments
    NM = sum(NProgMoms)
    mom_norms::NTuple{NM, FT} = CL.get_moments_normalizing_factors(Int.(NProgMoms), norms)
    
    # Define collision kernel function
    kernel_func = CL.KernelFunctions.LinearKernelFunction(5e0) # 5 m^3 / kg / s
    # kernel_func = CL.KernelFunctions.LongKernelFunction(5.236e-10, 9.44e9, 5.78)

    # Compute matrix of kernel tensors each approximating kernel func as polynomail serries
    r = 2
    matrix_of_kernels = ntuple(ND) do i
        ntuple(ND) do j
            if i == j == 1
                CL.KernelTensors.CoalescenceTensor(kernel_func, r, FT(5e-10))
            else
                CL.KernelTensors.CoalescenceTensor(kernel_func, r, FT(1e-6), FT(5e-10))
            end
        end
    end

    # Define mass thresholds between water categories
    mass_thresholds = (5e-10, Inf) # 5e-10 kg

    # Define coalescence data required by Cloudy
    coal_data::CL.Coalescence.CoalescenceData{ND, r+1, FT, (r+1)^2} = CL.Coalescence.CoalescenceData(matrix_of_kernels, NProgMoms, mass_thresholds, norms)
    
    # Define terminal velocity coefficients, assuming vt = sum_i v_i[1] * x^(v_i[2]) 
    # v1 is normalized by mass norm; v1 = v1 * norm[2] ^ v2
    vel = ((FT(30), FT(1.0 / 6)),) # 30 kg ^ (-1/6) * m / s
    cloudy_params = CO.Parameters.CloudyParameters(NProgMoms, norms, mom_norms, coal_data, vel)
    return cloudy_params, pdists
end
#! format: on
