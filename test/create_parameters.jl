import Kinematic1D as KID
import CLIMAParameters as CP
import CloudMicrophysics as CM
import Thermodynamics as TD

#! format: off
function create_parameter_set(
    out_dir::String,
    toml_dict::CP.AbstractTOMLDict,
    FTD = CP.float_type(toml_dict),
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
    z_0 = 0.0,
    z_1 = 740.0,
    z_2 = 3260.0,
    rv_0 = 0.015,
    rv_1 = 0.0138,
    rv_2 = 0.0024,
    tht_0 = 297.9,
    tht_1 = 297.9,
    tht_2 = 312.66,
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
        println(io, "[init_cond_z0]")
        println(io, "alias = \"z_0\"")
        println(io, "value = " * string(z_0))
        println(io, "type = \"float\"")
        println(io, "[init_cond_z1]")
        println(io, "alias = \"z_1\"")
        println(io, "value = " * string(z_1))
        println(io, "type = \"float\"")
        println(io, "[init_cond_z2]")
        println(io, "alias = \"z_2\"")
        println(io, "value = " * string(z_2))
        println(io, "type = \"float\"")
        println(io, "[init_cond_rv0]")
        println(io, "alias = \"rv_0\"")
        println(io, "value = " * string(rv_0))
        println(io, "type = \"float\"")
        println(io, "[init_cond_rv1]")
        println(io, "alias = \"rv_1\"")
        println(io, "value = " * string(rv_1))
        println(io, "type = \"float\"")
        println(io, "[init_cond_rv2]")
        println(io, "alias = \"rv_2\"")
        println(io, "value = " * string(rv_2))
        println(io, "type = \"float\"")
        println(io, "[init_cond_theta0]")
        println(io, "alias = \"tht_0\"")
        println(io, "value = " * string(tht_0))
        println(io, "type = \"float\"")
        println(io, "[init_cond_theta1]")
        println(io, "alias = \"tht_1\"")
        println(io, "value = " * string(tht_1))
        println(io, "type = \"float\"")
        println(io, "[init_cond_theta2]")
        println(io, "alias = \"tht_2\"")
        println(io, "value = " * string(tht_2))
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

    aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TD.Parameters.ThermodynamicsParameters{FTD}(; param_pairs...)
    TP = typeof(thermo_params)

    aliases = string.(fieldnames(CM.Parameters.ModalNucleationParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    modal_nucleation_params = CM.Parameters.ModalNucleationParameters{FT}(; pairs...)
    MNP = typeof(modal_nucleation_params)

    aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
    aliases = setdiff(aliases, ["thermo_params, modal_nucleation_params"])
    pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    microphys_params = CM.Parameters.CloudMicrophysicsParameters{FTD, TP, MNP}(;
        pairs...,
        thermo_params,
        modal_nucleation_params,
    )
    MP = typeof(microphys_params)

    aliases = ["w1", "t1", "p0", "precip_sources", "precip_sinks", "qtot_flux_correction", "prescribed_Nd", "r_dry", "std_dry", "κ", "z_0", "z_1", "z_2", "rv_0", "rv_1", "rv_2", "tht_0", "tht_1", "tht_2"]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Kinematic1D")

    param_set = KID.Parameters.KinematicParameters{FTD, MP}(; pairs..., microphys_params)
    if !isbits(param_set)
        print(param_set)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return param_set
end
#! format: on
