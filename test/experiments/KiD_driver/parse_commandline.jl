"""
    Parse the command line arguments for KiD_driver.jl
"""

function parse_commandline()
    s = AP.ArgParseSettings()

    AP.@add_arg_table! s begin
        "--moisture_choice"
        help = "Mositure model choice: EquilibriumMoisture, NonEquilibriumMoisture"
        arg_type = String
        default = "NonEquilibriumMoisture"
        "--prognostic_vars"
        help = "Prognostic variables choice: RhoThetaQ, RhodTQ"
        arg_type = String
        default = "RhoThetaQ"
        "--precipitation_choice"
        help = "Precipitation model choice: NoPrecipitation, Precipitation0M, Precipitation1M, Precipitation2M"
        arg_type = String
        default = "Precipitation1M"
        "--rain_formation_scheme_choice"
        help = "Rain formation scheme choice: CliMA_1M, KK2000, B1994, TC1980, LD2004 for Precipitation1M; and SB2006 for Precipitation2M"
        arg_type = String
        default = "CliMA_1M"
        "--prescribed_Nd"
        help = "Prescribed number of cloud droplets (used in KK2000, B1994, TC1980, LD2004 and SB2006 rain formation schemes)"
        arg_type = Float64
        default = Float64(1e8)
        "--plotting_flag"
        help = "Set to true if you want to generate some basic plots at the end of the simulation"
        arg_type = Bool
        default = true
        "--precip_sources"
        help = "Set to true if you want to switch on autoconversion and accretion in the 1-moment scheme"
        arg_type = Bool
        default = true
        "--precip_sinks"
        help = "Set to true if you want to switch on evaporation, deposition, sublimation and melting in the 1-moment scheme"
        arg_type = Bool
        default = true
        "--qtot_flux_correction"
        help = "Set to true if you want to apply flux correction for advecting q_tot.
        (By default flux correction is not applied to q_tot but is applied to all other microphysics tracers)"
        arg_type = Bool
        default = false
        "--z_min"
        help = "Bottom of the computational domain [m]"
        arg_type = Real
        default = 0.0
        "--z_max"
        help = "Top of the computational domain [m]"
        arg_type = Real
        default = 2000.0
        "--n_elem"
        help = "Number of computational elements"
        arg_type = Int
        default = 60
        "--dt"
        help = "Simulation time step [s]"
        arg_type = Real
        default = 1.0
        "--dt_output"
        help = "Output time step [s]"
        arg_type = Real
        default = 30.0
        "--t_ini"
        help = "Time at the beginning of the simulation [s]"
        arg_type = Real
        default = 0.0
        "--t_end"
        help = "Time at the end of the simulation [s]"
        arg_type = Real
        default = 3600.0
        "--w1"
        help = "Maximum prescribed updraft momentum flux [m/s * kg/m3]"
        arg_type = Real
        default = 2.0
        "--t1"
        help = "Oscillation time of the prescribed momentum flux [s]"
        arg_type = Real
        default = 600.0
        "--p0"
        help = "Pressure at the surface [pa]"
        arg_type = Real
        default = 99400.0 #100000.0
        "--r_dry"
        help = "aerosol distribution mean radius for aerosol activation calculations in 2M schemes [m]"
        arg_type = Real
        default = 0.04 * 1e-6
        "--std_dry"
        help = "aerosol distribution standard deviation for aerosol activation calucaulations in 2M schemes"
        arg_type = Real
        default = 1.4
        "--kappa"
        help = "hygroscopicity of aerosols for aerosol activation calucaulations in 2M schemes"
        arg_type = Real
        default = 0.9
    end

    return AP.parse_args(s)
end
