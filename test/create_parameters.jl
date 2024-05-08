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
    κ = 0.9,
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
        "prescribed_flow_w1" => Dict("value" => w1, "type" => "float"),
        "prescribed_flow_t1" => Dict("value" => t1, "type" => "float"),
        "surface_pressure" => Dict("value" => p0, "type" => "float"),
        "precipitation_sources_flag" => Dict("value" => precip_sources, "type" => "bool"),
        "precipitation_sinks_flag" => Dict("value" => precip_sinks, "type" => "bool"),
        "qtot_flux_correction_flag" => Dict("value" => qtot_flux_correction, "type" => "bool"),
        "prescribed_Nd" => Dict("value" => prescribed_Nd, "type" => "float"),
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
Base.broadcastable(x::CO.Parameters.CommonParameters) = Ref(x)

function create_kid_parameters(toml_dict)
    kid_params = KID.Parameters.KinematicDriverParameters(toml_dict)
    if !isbits(kid_params)
        print(kid_params)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return kid_params
end
Base.broadcastable(x::KID.Parameters.KinematicDriverParameters) = Ref(x)

function create_cloudy_parameters(NM)
    pdists::NTuple{Int(NM / 3), CL.ParticleDistributions.GammaPrimitiveParticleDistribution{FT}} = ntuple(Int(NM/3)) do k
        CL.ParticleDistributions.GammaPrimitiveParticleDistribution(FT(0), FT(k), FT(1))
    end # TODO: generalize
    NProgMoms::NTuple{Int(NM / 3), FT} = ntuple(length(pdists)) do k
        CL.ParticleDistributions.nparams(pdists[k])
    end
    # TODO: figure out why norms isn't working
    norms::Tuple{FT, FT} = (FT(1), FT(1)) #(FT(1e6), FT(1e-9))
    mom_norms::NTuple{NM, FT} = CL.get_moments_normalizing_factors(Int.(NProgMoms), norms)
    
    kernel_func = CL.KernelFunctions.LongKernelFunction(5.236e-10, 9.44e9, 5.78)
    matrix_of_kernels = ntuple(2) do i
        ntuple(2) do j
            if i == j == 1
                CL.KernelTensors.CoalescenceTensor(kernel_func, 2, FT(5e-10))
            else
                CL.KernelTensors.CoalescenceTensor(kernel_func, 2, FT(1e-6), FT(5e-10))
            end
        end
    end
    coal_data::CL.Coalescence.CoalescenceData{2, 3, FT, 9} = CL.Coalescence.CoalescenceData(matrix_of_kernels, Int.(NProgMoms), (5e-10, Inf), norms)
    vel = ((FT(50.0), FT(1.0 / 6)),)
    cloudy_params = CO.Parameters.CloudyParameters(NProgMoms, norms, mom_norms, coal_data, vel)
    return cloudy_params
end
Base.broadcastable(x::CO.Parameters.CloudyParameters) = Ref(x)
#! format: on
