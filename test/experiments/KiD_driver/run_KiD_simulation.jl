import NCDatasets as NC
import OrdinaryDiffEq as ODE
import UnPack
import ArgParse as AP

import ClimaCore as CC
import Thermodynamics as TD
import CloudMicrophysics.CommonTypes as CMT
import CLIMAParameters as CP
import Kinematic1D as KID

const kid_dir = pkgdir(KID)
include(joinpath(kid_dir, "test", "create_parameters.jl"))

include("../../plotting_utils.jl")


function run_KiD_simulation(::Type{FT}, opts) where {FT}

    # Equations to solve for mositure and precipitation variables
    moisture_choice = opts["moisture_choice"]
    prognostics_choice = opts["prognostic_vars"]
    precipitation_choice = opts["precipitation_choice"]
    rain_formation_choice = opts["rain_formation_scheme_choice"]

    if moisture_choice == "EquilibriumMoisture" && prognostics_choice == "RhoThetaQ"
        moisture = KID.EquilibriumMoisture_ρθq()
    elseif moisture_choice == "EquilibriumMoisture" && prognostics_choice == "RhodTQ"
        moisture = KID.EquilibriumMoisture_ρdTq()
    elseif moisture_choice == "NonEquilibriumMoisture" && prognostics_choice == "RhoThetaQ"
        moisture = KID.NonEquilibriumMoisture_ρθq()
    elseif moisture_choice == "NonEquilibriumMoisture" && prognostics_choice == "RhodTQ"
        moisture = KID.NonEquilibriumMoisture_ρdTq()
    else
        error("Invalid moisture choice: $moisture_choice or invalid prognostic variables choice: $prognostic_vars")
    end
    if precipitation_choice == "NoPrecipitation"
        precip = KID.NoPrecipitation()
    elseif precipitation_choice == "Precipitation0M"
        precip = KID.Precipitation0M()
    elseif precipitation_choice == "Precipitation1M"
        if rain_formation_choice == "CliMA_1M"
            precip = KID.Precipitation1M(KID.OneMomentRainFormation())
        elseif rain_formation_choice == "KK2000"
            precip = KID.Precipitation1M(CMT.KK2000Type())
        elseif rain_formation_choice == "B1994"
            precip = KID.Precipitation1M(CMT.B1994Type())
        elseif rain_formation_choice == "TC1980"
            precip = KID.Precipitation1M(CMT.TC1980Type())
        elseif rain_formation_choice == "LD2004"
            precip = KID.Precipitation1M(CMT.LD2004Type())
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
    elseif precipitation_choice == "Precipitation2M"
        if rain_formation_choice == "SB2006"
            precip = KID.Precipitation2M(CMT.SB2006Type())
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
    else
        error("Invalid precipitation choice: $precipitation_choice")
    end

    # Initialize the timestepping struct
    TS = KID.TimeStepping(FT(opts["dt"]), FT(opts["dt_output"]), FT(opts["t_end"]))

    # Create the coordinates
    space, face_space = KID.make_function_space(FT, FT(opts["z_min"]), FT(opts["z_max"]), opts["n_elem"])
    coord = CC.Fields.coordinate_field(space)
    face_coord = CC.Fields.coordinate_field(face_space)

    # Initialize the netcdf output Stats struct
    output_folder = string("Output_", moisture_choice, "_", prognostics_choice, "_", precipitation_choice)
    if precipitation_choice in ["Precipitation1M", "Precipitation2M"]
        output_folder = output_folder * "_" * rain_formation_choice
    end
    if opts["qtot_flux_correction"]
        output_folder = output_folder * "_wFC"
    end
    path = joinpath(@__DIR__, output_folder)
    mkpath(path)
    fname = joinpath(path, "Output.nc")
    Stats = KID.NetCDFIO_Stats(fname, 1.0, vec(face_coord), vec(coord))

    # Instantiate CliMA Parameters and overwrite the defaults
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    params = create_parameter_set(
        path,
        toml_dict,
        FT,
        FT(opts["w1"]),
        FT(opts["t1"]),
        FT(opts["p0"]),
        Int(opts["precip_sources"]),
        Int(opts["precip_sinks"]),
        Int(opts["qtot_flux_correction"]),
        FT(opts["prescribed_Nd"]),
        FT(opts["r_dry"]),
        FT(opts["std_dry"]),
        FT(opts["kappa"]),
    )

    # Solve the initial value problem for density profile
    ρ_profile = KID.ρ_ivp(FT, params)
    # Create the initial condition profiles
    init = map(coord -> KID.init_1d_column(FT, params, ρ_profile, coord.z), coord)

    # Create state vector and apply initial condition
    Y = KID.initialise_state(moisture, precip, init)

    # Create aux vector and apply initial condition
    aux = KID.initialise_aux(FT, init, params, TS, Stats, face_space, moisture)

    # Output the initial condition
    KID.KiD_output(aux, 0.0)

    # Define callbacks for output
    callback_io = ODE.DiscreteCallback(KID.condition_io, KID.affect_io!; save_positions = (false, false))
    callbacks = ODE.CallbackSet(callback_io)

    # Collect all the tendencies into rhs function for ODE solver
    # based on model choices for the solved equations
    ode_rhs! = KID.make_rhs_function(moisture, precip)

    # Solve the ODE operator
    problem = ODE.ODEProblem(ode_rhs!, Y, (FT(opts["t_ini"]), FT(opts["t_end"])), aux)
    solver = ODE.solve(
        problem,
        ODE.SSPRK33(),
        dt = TS.dt,
        saveat = TS.dt_io,
        callback = callbacks,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
    )

    # Some basic plots
    if opts["plotting_flag"] == true

        plot_folder = string("experiments/KiD_driver/", output_folder, "/figures/")

        z_centers = parent(CC.Fields.coordinate_field(space))
        plot_final_aux_profiles(z_centers, aux, precip, output = plot_folder)
        plot_animation(z_centers, solver, aux, moisture, precip, output = plot_folder)
        plot_timeheight(string("experiments/KiD_driver/", output_folder, "/Output.nc"), output = plot_folder)
    end

    return solver

end
