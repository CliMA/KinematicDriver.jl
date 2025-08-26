import KinematicDriver.Common as CO
import KinematicDriver.K1DModel as KID
import ClimaParams as CP
import CloudMicrophysics as CM
import Thermodynamics as TD
import Cloudy as CL

function override_toml_dict(
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
    open_system_activation = 0,
    local_activation = 0,
    r_dry = 0.04 * 1e-6,
    std_dry = 1.4,
    κ = 1.12,
)
    FT = CP.float_type(toml_dict)
    override_file = Dict(
        "mean_sea_level_pressure" => Dict("value" => 100000.0, "type" => "float"),
        "gravitational_acceleration" => Dict("value" => 9.80665, "type" => "float"),
        "universal_gas_constant" => Dict("value" => 8.314462618, "type" => "float"),
        "adiabatic_exponent_dry_air" => Dict("value" => 0.2855747338575384, "type" => "float"),
        "isobaric_specific_heat_vapor" => Dict("value" => 1850.0, "type" => "float"),
        "gas_constant_vapor" => Dict("value" => 461.52998157091315, "type" => "float"),
        "gas_constant_dry_air" => Dict("value" => 287.0027047999343, "type" => "float"),
        "cloud_liquid_water_specific_humidity_autoconversion_threshold" =>
            Dict("value" => 0.0001, "type" => "float"),
        "prescribed_flow_w1" => Dict("value" => w1, "type" => "float"),
        "prescribed_flow_t1" => Dict("value" => t1, "type" => "float"),
        "surface_pressure" => Dict("value" => p0, "type" => "float"),
        "precipitation_sources_flag" => Dict("value" => precip_sources, "type" => "bool"),
        "precipitation_sinks_flag" => Dict("value" => precip_sinks, "type" => "bool"),
        "qtot_flux_correction_flag" => Dict("value" => qtot_flux_correction, "type" => "bool"),
        "prescribed_Nd" => Dict("value" => prescribed_Nd, "type" => "float"),
        "open_system_activation" => Dict("value" => open_system_activation, "type" => "bool"),
        "local_activation" => Dict("value" => local_activation, "type" => "bool"),
        "r_dry" => Dict("value" => r_dry, "type" => "float"),
        "std_dry" => Dict("value" => std_dry, "type" => "float"),
        "kappa" => Dict("value" => κ, "type" => "float"),
        "init_cond_z0" => Dict("value" => z_0, "type" => "float"),
        "init_cond_z1" => Dict("value" => z_1, "type" => "float"),
        "init_cond_z2" => Dict("value" => z_2, "type" => "float"),
        "init_cond_rv0" => Dict("value" => rv_0, "type" => "float"),
        "init_cond_rv1" => Dict("value" => rv_1, "type" => "float"),
        "init_cond_rv2" => Dict("value" => rv_2, "type" => "float"),
        "init_cond_theta0" => Dict("value" => tht_0, "type" => "float"),
        "init_cond_theta1" => Dict("value" => tht_1, "type" => "float"),
        "init_cond_theta2" => Dict("value" => tht_2, "type" => "float"),
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

function determine_cloudy_disttypes(NM::Int, default_gamma = true, initial_gamma = false)
    if NM < 2
        error("At least two moments are required")
    end
    if default_gamma
        if NM % 3 == 0
            ND = Int(NM / 3)
            dist_names = ntuple(_ -> "gamma", ND)
        elseif NM % 3 == 2
            ND = Int((NM - 2) / 3 + 1)
            expmode = initial_gamma ? ND : 1
            dist_names = ntuple(ND) do k
                if k == expmode
                    "exponential"
                else
                    "gamma"
                end
            end
        else
            ND = Int((NM - 4) / 3 + 2)
            expmodes = initial_gamma ? (ND - 1, ND) : (1, 2)
            dist_names = ntuple(ND) do k
                if k <= expmodes[2] && k >= expmodes[1]
                    "exponential"
                else
                    "gamma"
                end
            end
        end
    else # default exponential
        if NM % 2 == 0
            ND = Int(NM / 2)
            dist_names = ntuple(_ -> "exponential", ND)
        else
            ND = Int((NM - 3) / 2 + 1)
            gmode = initial_gamma ? 1 : ND
            dist_names = ntuple(ND) do k
                if k == gmode
                    "gamma"
                else
                    "exponential"
                end
            end
        end
    end
    @show dist_names
    return dist_names
end

function create_cloudy_parameters(FT, dist_names::NTuple{ND, String} = ("gamma", "gamma")) where {ND}

    # Create water category distributions
    function make_dist(dist_name)
        if dist_name == "monodisperse"
            return CL.ParticleDistributions.MonodispersePrimitiveParticleDistribution(FT(0), FT(1))
        elseif dist_name == "exponential"
            return CL.ParticleDistributions.ExponentialPrimitiveParticleDistribution(FT(0), FT(1))
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
    size_threshold = FT(5e-10)
    mass_thresholds = ntuple(ND) do k
        if k < ND
            size_threshold * 10^(k - ND / 2)
        else
            Inf
        end
    end
    @show mass_thresholds

    # Define coalescence data required by Cloudy
    coal_data::CL.Coalescence.CoalescenceData{ND, r + 1, FT, (r + 1)^2} =
        CL.Coalescence.CoalescenceData(matrix_of_kernels, NProgMoms, mass_thresholds, norms)

    # Define terminal velocity coefficients, assuming vt = sum_i v_i[1] * x^(v_i[2]) 
    # v1 is normalized by mass norm; v1 = v1 * norm[2] ^ v2
    vel = ((FT(30), FT(1.0 / 6)),) # 30 kg ^ (-1/6) * m / s

    cloudy_params = CO.Parameters.CloudyParameters(NProgMoms, norms, mom_norms, coal_data, vel, size_threshold)
    return cloudy_params, pdists
end
