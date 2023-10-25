import NCDatasets as NC
import OrdinaryDiffEq as ODE

import ClimaCore as CC
import CLIMAParameters as CP
import CloudMicrophysics.Parameters as CMP
import Kinematic1D as KID

include(joinpath(pkgdir(KID), "test", "create_parameters.jl"))
include(joinpath(pkgdir(KID), "test", "plotting_utils.jl"))

function run_KiD_simulation(::Type{FT}, opts) where {FT}

    # Equations to solve for mositure and precipitation variables
    moisture_choice = opts["moisture_choice"]
    prognostics_choice = opts["prognostic_vars"]
    precipitation_choice = opts["precipitation_choice"]
    rain_formation_choice = opts["rain_formation_scheme_choice"]
    sedimentation_choice = opts["sedimentation_scheme_choice"]

    # Decide the output flder name based on options
    output_folder = string("Output_", moisture_choice, "_", prognostics_choice, "_", precipitation_choice)
    if precipitation_choice in ["Precipitation1M", "Precipitation2M"]
        output_folder = output_folder * "_" * rain_formation_choice
        if sedimentation_choice == "Chen2022"
            output_folder = output_folder * "_Chen2022"
        end
    end
    if opts["qtot_flux_correction"]
        output_folder = output_folder * "_wFC"
    end
    path = joinpath(@__DIR__, output_folder)
    mkpath(path)

    # Overwrite the defaults parameters based on options
    default_toml_dict = CP.create_toml_dict(FT, dict_type = "alias")
    toml_dict = override_toml_dict(
        path,
        default_toml_dict,
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
    # Create Thermodynamics.jl and Kinematic1D model parameters
    # (some of the CloudMicrophysics.jl parameters structs are created later based on model choices)
    thermo_params = create_thermodynamics_parameters(toml_dict)
    kid_params = create_kid_parameters(toml_dict)
    air_params = CMP.AirProperties(FT, toml_dict)
    activation_params = CMP.AerosolActivationParameters(FT, toml_dict)

    if moisture_choice == "EquilibriumMoisture" && prognostics_choice == "RhoThetaQ"
        moisture = KID.EquilibriumMoisture_ρθq()
    elseif moisture_choice == "EquilibriumMoisture" && prognostics_choice == "RhodTQ"
        moisture = KID.EquilibriumMoisture_ρdTq()
    elseif moisture_choice == "NonEquilibriumMoisture" && prognostics_choice == "RhoThetaQ"
        moisture = KID.NonEquilibriumMoisture_ρθq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
    elseif moisture_choice == "NonEquilibriumMoisture" && prognostics_choice == "RhodTQ"
        moisture = KID.NonEquilibriumMoisture_ρdTq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
    else
        error("Invalid moisture choice: $moisture_choice or invalid prognostic variables choice: $prognostic_vars")
    end

    if precipitation_choice == "NoPrecipitation"
        precip = KID.NoPrecipitation()
    elseif precipitation_choice == "Precipitation0M"
        precip = KID.Precipitation0M(CMP.Parameters0M(FT, toml_dict))
    elseif precipitation_choice == "Precipitation1M"
        if sedimentation_choice == "CliMA_1M"
            st = CMP.Blk1MVelType(FT, toml_dict)
        elseif sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelType(FT, toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "CliMA_1M"
            rain_params = CMP.Rain(FT, toml_dict)
            rf = rain_params.acnv1M
        elseif rain_formation_choice == "KK2000"
            rf = CMP.KK2000(FT, toml_dict)
        elseif rain_formation_choice == "B1994"
            rf = CMP.B1994(FT, toml_dict)
        elseif rain_formation_choice == "TC1980"
            rf = CMP.TC1980(FT, toml_dict)
        elseif rain_formation_choice == "LD2004"
            rf = CMP.LD2004(FT, toml_dict)
        elseif rain_formation_choice == "VarTimeScaleAcnv"
            rf = CMP.VarTimescaleAcnv(FT, toml_dict)
        else
            error("Invalid rain formation choice: $rain_formation_choice")
        end
        precip = KID.Precipitation1M(
            CMP.CloudLiquid(FT, toml_dict),
            CMP.CloudIce(FT, toml_dict),
            CMP.Rain(FT, toml_dict),
            CMP.Snow(FT, toml_dict),
            CMP.CollisionEff(FT, toml_dict),
            rf,
            st,
        )
    elseif precipitation_choice == "Precipitation2M"
        if sedimentation_choice == "Chen2022"
            st = CMP.Chen2022VelTypeRain(FT, toml_dict)
        elseif sedimentation_choice == "SB2006"
            st = CMP.SB2006VelType(FT, toml_dict)
        else
            error("Invalid sedimentation choice: $sedimentation_choice")
        end
        if rain_formation_choice == "SB2006"
            precip = KID.Precipitation2M(CMP.SB2006(FT, toml_dict), st)
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
    fname = joinpath(path, "Output.nc")
    Stats = KID.NetCDFIO_Stats(fname, 1.0, parent(face_coord), parent(coord))

    # Solve the initial value problem for density profile
    ρ_profile = KID.ρ_ivp(FT, kid_params, thermo_params)
    # Create the initial condition profiles
    init = map(coord -> KID.init_1d_column(FT, kid_params, thermo_params, ρ_profile, coord.z), coord)

    # Create state vector and apply initial condition
    Y = KID.initialise_state(moisture, precip, init)

    # Create aux vector and apply initial condition
    aux = KID.initialise_aux(
        FT,
        init,
        kid_params,
        thermo_params,
        air_params,
        activation_params,
        TS,
        Stats,
        face_space,
        moisture,
    )

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
        plot_animation(z_centers, solver, aux, moisture, precip, KID, output = plot_folder)
        plot_timeheight(string("experiments/KiD_driver/", output_folder, "/Output.nc"), output = plot_folder)
    end

    return solver

end
