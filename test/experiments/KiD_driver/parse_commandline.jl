"""
    Parse the command line arguments for KiD_driver.jl
"""

import ArgParse as AP

function parse_commandline()
    s = AP.ArgParseSettings()

    AP.@add_arg_table! s begin
        "--FLOAT_TYPE"
        help = "Float type. Can be set to Float64 or Float32"
        arg_type = String
        default = "Float64"
        "--moisture_choice"
        help = "Mositure model choice: EquilibriumMoisture, NonEquilibriumMoisture, CloudyMoisture, MoistureP3"
        arg_type = String
        default = "EquilibriumMoisture"
        "--precipitation_choice"
        help = "Precipitation model choice: NoPrecipitation, Precipitation0M, Precipitation1M, Precipitation2M, CloudyPrecip, PrecipitationP3"
        arg_type = String
        default = "Precipitation1M"
        "--num_moments"
        help = "Number of moments to use for CloudyPrecip (ignored otherwise)"
        arg_type = Int
        default = 6
        "--rain_formation_scheme_choice"
        help = "Rain formation scheme choice: CliMA_1M, KK2000, B1994, TC1980, LD2004, VarTimeScaleAcnv for Precipitation1M; and SB2006, SB2006NL for Precipitation2M"
        arg_type = String
        default = "CliMA_1M"
        "--sedimentation_scheme_choice"
        help = "Sedimentation scheme choice: CliMA_1M, Chen2022 for Precipitation1M; and Chen2022, SB2006 for Precipitation2M"
        arg_type = String
        default = "CliMA_1M"
        "--prescribed_Nd"
        help = "Prescribed number of cloud droplets (used in KK2000, B1994, TC1980, LD2004, VarTimeScaleAcnv and SB2006 rain formation schemes)"
        arg_type = Float64
        default = Float64(1e8)
        "--open_system_activation"
        help = "Set to true if you want to assume an open system for aerosol activation with an aerosol concentration budget equal to prescribed_Nd"
        arg_type = Bool
        default = false
        "--local_activation"
        help = "Set to true if you want to apply ARG aerosol activation locally otherwise activation occurs at the cloud base"
        arg_type = Bool
        default = false
        "--plotting_flag"
        help = "Set to true if you want to generate some basic plots at the end of the simulation"
        arg_type = Bool
        default = true
        "--precip_sources"
        help = "Set to true if you want to switch on autoconversion and accretion in the 1-moment scheme, or collisional coalescence in Cloudy"
        arg_type = Bool
        default = true
        "--precip_sinks"
        help = "Set to true if you want to switch on evaporation, deposition, sublimation and melting in the 1-moment scheme; or condensation/evaporation in Cloudy"
        arg_type = Bool
        default = true
        "--qtot_flux_correction"
        help = "Set to true if you want to apply flux correction for advecting q_tot.
        (By default flux correction is not applied to q_tot but is applied to all other microphysics tracers)"
        arg_type = Bool
        default = false
        "--z_min"
        help = "Bottom of the computational domain [m]"
        arg_type = Float64
        default = Float64(0.0)
        "--z_max"
        help = "Top of the computational domain [m]"
        arg_type = Float64
        default = Float64(2000)
        "--n_elem"
        help = "Number of computational elements"
        arg_type = Int
        default = 128
        "--dt"
        help = "Simulation time step [s]"
        arg_type = Float64
        default = Float64(1)
        "--dt_output"
        help = "Output time step [s]"
        arg_type = Float64
        default = Float64(30)
        "--t_ini"
        help = "Time at the beginning of the simulation [s]"
        arg_type = Float64
        default = Float64(0)
        "--t_end"
        help = "Time at the end of the simulation [s]"
        arg_type = Float64
        default = Float64(3600)
        "--w1"
        help = "Maximum prescribed updraft momentum flux [m/s * kg/m3]"
        arg_type = Float64
        default = Float64(2)
        "--t1"
        help = "Oscillation time of the prescribed momentum flux [s]"
        arg_type = Float64
        default = Float64(600)
        "--p0"
        help = "Pressure at the surface [pa]"
        arg_type = Float64
        default = Float64(100000)
        "--r_dry"
        help = "aerosol distribution mean radius for aerosol activation calculations in 2M schemes [m]"
        arg_type = Float64
        default = Float64(0.04 * 1e-6)
        "--std_dry"
        help = "aerosol distribution standard deviation for aerosol activation calucaulations in 2M schemes"
        arg_type = Float64
        default = Float64(1.4)
        "--kappa"
        help = "hygroscopicity of aerosols for aerosol activation calucaulations in 2M schemes"
        arg_type = Float64
        default = Float64(0.9)
        "--z_0"
        help = "Initial condition z0 [m]"
        arg_type = Float64
        default = Float64(0)
        "--z_1"
        help = "Initial condition z1 [m]"
        arg_type = Float64
        default = Float64(740)
        "--z_2"
        help = "Initial condition z2 [m]"
        arg_type = Float64
        default = Float64(3260)
        "--rv_0"
        help = "Initial condition rv0 [kg/kg]"
        arg_type = Float64
        default = Float64(0.015)
        "--rv_1"
        help = "Initial condition rv1 [kg/kg]"
        arg_type = Float64
        default = Float64(0.0138)
        "--rv_2"
        help = "Initial condition rv2 [kg/kg]"
        arg_type = Float64
        default = Float64(0.0024)
        "--tht_0"
        help = "Initial condition theta0 [K]"
        arg_type = Float64
        default = Float64(297.9)
        "--tht_1"
        help = "Initial condition theta1 [K]"
        arg_type = Float64
        default = Float64(297.9)
        "--tht_2"
        help = "Initial condition theta2 [K]"
        arg_type = Float64
        default = Float64(312.66)
        "--p3_boundary_condition"
        help = "Characteristics of particle flux being introduced into P3 domain top (if ice_start = false) or of the initial ice signal (if ice_start = true)"
        arg_type = NamedTuple
        default = (;
            ice_start = false,
            _magnitude = Float64(0.5),
            _q_flux = Float64(0.65e-4),
            _N_flux = Float64(40000),
            _F_rim = Float64(0.2),
            _F_liq = Float64(0.2),
            _ρ_r_init = Float64(900),
        )
        "--init_sounding"
        help = "Exterior data used for KiD initial conditions"
        arg_type = String
        default = "ShipwayHill"

    end

    return AP.parse_args(s)
end
