import Kinematic1D as KID
import CLIMAParameters as CP
import CloudMicrophysics as CM
import Thermodynamics as TD

#! format: off
function override_toml_dict(
    out_dir::String,
    toml_dict::CP.AbstractTOMLDict,
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
        println(io, "value = 100000.0")
        println(io, "type = \"float\"")
        println(io, "[gravitational_acceleration]")
        println(io, "value = 9.80665")
        println(io, "type = \"float\"")
        println(io, "[gas_constant]")
        println(io, "value = 8.314462618")
        println(io, "type = \"float\"")
        println(io, "[adiabatic_exponent_dry_air]")
        println(io, "value = 0.2855747338575384")
        println(io, "type = \"float\"")
        println(io, "[isobaric_specific_heat_vapor]")
        println(io, "value = 1850.0")
        println(io, "type = \"float\"")
        println(io, "[molar_mass_dry_air]")
        println(io, "value = 0.02896998")
        println(io, "type = \"float\"")
        println(io, "[molar_mass_water]")
        println(io, "value = 0.018015")
        println(io, "type = \"float\"")
        println(io, "[cloud_liquid_water_specific_humidity_autoconversion_threshold]")
        println(io, "value = 0.0001")
        println(io, "type = \"float\"")
        println(io, "[prescribed_flow_w1]")
        println(io, "value = " * string(w1))
        println(io, "type = \"float\"")
        println(io, "[prescribed_flow_t1]")
        println(io, "value = " * string(t1))
        println(io, "type = \"float\"")
        println(io, "[surface_pressure]")
        println(io, "value = " * string(p0))
        println(io, "type = \"float\"")
        println(io, "[precipitation_sources_flag]")
        println(io, "value = " * string(precip_sources))
        println(io, "type = \"integer\"")
        println(io, "[precipitation_sinks_flag]")
        println(io, "value = " * string(precip_sinks))
        println(io, "type = \"integer\"")
        println(io, "[qtot_flux_correction_flag]")
        println(io, "value = " * string(qtot_flux_correction))
        println(io, "type = \"integer\"")
        println(io, "[prescribed_Nd]")
        println(io, "value = " * string(prescribed_Nd))
        println(io, "type = \"float\"")
        println(io, "[r_dry]")
        println(io, "value = " * string(r_dry))
        println(io, "type = \"float\"")
        println(io, "[std_dry]")
        println(io, "value = " * string(std_dry))
        println(io, "type = \"float\"")
        println(io, "[kappa]")
        println(io, "value = " * string(κ))
        println(io, "type = \"float\"")
    end
    toml_dict = CP.create_toml_dict(FT; override_file)
    isfile(override_file) && rm(override_file; force=true)
    return toml_dict
end

function create_thermodynamics_parameters(toml_dict)
    FTD = CP.float_type(toml_dict)
    thermo_params = TD.Parameters.ThermodynamicsParameters(toml_dict)
    return thermo_params
end

function create_kid_parameters(toml_dict)
    FTD = CP.float_type(toml_dict)
    kid_params = KID.Parameters.KinematicParameters(FTD, toml_dict)
    if !isbits(kid_params)
        print(kid_params)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return kid_params
end
Base.broadcastable(x::KID.Parameters.KinematicParameters) = Ref(x)
#! format: on
