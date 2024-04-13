import Kinematic1D.Common as CO
import Kinematic1D.K1DModel as KID
import CLIMAParameters as CP
import CloudMicrophysics as CM
import Thermodynamics as TD

#! format: off
function override_toml_dict(
    out_dir::String,
    toml_dict::CP.AbstractTOMLDict;
    w1 = 2.0,
    t1 = 600.0,
    p0 = 100700.0,
    precip_sources = 1,
    precip_sinks = 1,
    qtot_flux_correction = 0,
    prescribed_Nd = 100 * 1e6,
    r_dry = 0.04 * 1e-6,
    std_dry = 1.4,
    κ = 0.9,
)
    FT = CP.float_type(toml_dict)
    override_file = joinpath(out_dir, "override_dict.toml")
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
        println(io, "value = 0.0001")
        println(io, "type = \"float\"")
        println(io, "[prescribed_flow_w1]")
        println(io, "alias = \"w1\"")
        println(io, "value = " * string(w1))
        println(io, "type = \"float\"")
        println(io, "[prescribed_flow_t1]")
        println(io, "alias = \"t1\"")
        println(io, "value = " * string(t1))
        println(io, "type = \"float\"")
        println(io, "[surface_pressure]")
        println(io, "alias = \"p0\"")
        println(io, "value = " * string(p0))
        println(io, "type = \"float\"")
        println(io, "[precipitation_sources_flag]")
        println(io, "alias = \"precip_sources\"")
        println(io, "value = " * string(precip_sources))
        println(io, "type = \"integer\"")
        println(io, "[precipitation_sinks_flag]")
        println(io, "alias = \"precip_sinks\"")
        println(io, "value = " * string(precip_sinks))
        println(io, "type = \"integer\"")
        println(io, "[qtot_flux_correction_flag]")
        println(io, "alias = \"qtot_flux_correction\"")
        println(io, "value = " * string(qtot_flux_correction))
        println(io, "type = \"integer\"")
        println(io, "[prescribed_Nd]")
        println(io, "alias = \"prescribed_Nd\"")
        println(io, "value = " * string(prescribed_Nd))
        println(io, "type = \"float\"")
        println(io, "[r_dry]")
        println(io, "alias = \"r_dry\"")
        println(io, "value = " * string(r_dry))
        println(io, "type = \"float\"")
        println(io, "[std_dry]")
        println(io, "alias = \"std_dry\"")
        println(io, "value = " * string(std_dry))
        println(io, "type = \"float\"")
        println(io, "[kappa]")
        println(io, "alias = \"κ\"")
        println(io, "value = " * string(κ))
        println(io, "type = \"float\"")
    end
    toml_dict = CP.create_toml_dict(FT; override_file, dict_type="alias")
    isfile(override_file) && rm(override_file; force=true)
    return toml_dict
end

function create_thermodynamics_parameters(toml_dict)
    FTD = CP.float_type(toml_dict)
    aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TD.Parameters.ThermodynamicsParameters{FTD}(; param_pairs...)
    return thermo_params
end

function create_common_parameters(toml_dict)
    FTD = CP.float_type(toml_dict)
    aliases = ["precip_sources", "precip_sinks", "prescribed_Nd"]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Common")
    common_params = CO.Parameters.CommonParameters{FTD}(; pairs...)
    if !isbits(common_params)
        print(common_params)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return common_params
end
Base.broadcastable(x::CO.Parameters.CommonParameters) = Ref(x)

function create_kid_parameters(toml_dict)
    FTD = CP.float_type(toml_dict)
    aliases = ["w1", "t1", "p0", "qtot_flux_correction", "r_dry", "std_dry", "κ"]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Kinematic1D")
    kid_params = KID.Parameters.Kinematic1DParameters{FTD}(; pairs...)
    if !isbits(kid_params)
        print(kid_params)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return kid_params
end
Base.broadcastable(x::KID.Parameters.Kinematic1DParameters) = Ref(x)
#! format: on
