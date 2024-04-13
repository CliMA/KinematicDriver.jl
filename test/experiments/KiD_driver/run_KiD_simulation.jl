import OrdinaryDiffEq as ODE
import ClimaCore as CC
import CLIMAParameters as CP
import CloudMicrophysics.Parameters as CMP
import Kinematic1D
import Kinematic1D.Common as CO
import Kinematic1D.K1DModel as K1D

include(joinpath(pkgdir(Kinematic1D), "test", "create_parameters.jl"))
include(joinpath(pkgdir(Kinematic1D), "test", "plotting_utils.jl"))

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
        w1 = FT(opts["w1"]),
        t1 = FT(opts["t1"]),
        p0 = FT(opts["p0"]),
        precip_sources = Int(opts["precip_sources"]),
        precip_sinks = Int(opts["precip_sinks"]),
        qtot_flux_correction = Int(opts["qtot_flux_correction"]),
        prescribed_Nd = FT(opts["prescribed_Nd"]),
        r_dry = FT(opts["r_dry"]),
        std_dry = FT(opts["std_dry"]),
        kappa = FT(opts["kappa"]),
    )
    # Create Thermodynamics.jl and Kinematic1D model parameters
    # (some of the CloudMicrophysics.jl parameters structs are created later based on model choices)
    common_params = create_common_parameters(toml_dict)
    kid_params = create_kid_parameters(toml_dict)
    thermo_params = create_thermodynamics_parameters(toml_dict)
    air_params = CMP.AirProperties(FT, toml_dict)
    activation_params = CMP.AerosolActivationParameters(FT, toml_dict)

    moisture = CO.get_moisture_type(FT, moisture_choice, toml_dict)
    precip = CO.get_precipitation_type(
        FT,
        precipitation_choice,
        toml_dict;
        rain_formation_choice = rain_formation_choice,
        sedimentation_choice = sedimentation_choice,
    )

    # Initialize the timestepping struct
    TS = CO.TimeStepping(FT(opts["dt"]), FT(opts["dt_output"]), FT(opts["t_end"]))

    # Create the coordinates
    space, face_space = K1D.make_function_space(FT, FT(opts["z_min"]), FT(opts["z_max"]), opts["n_elem"])
    coord = CC.Fields.coordinate_field(space)
    face_coord = CC.Fields.coordinate_field(face_space)

    # Initialize the netcdf output Stats struct
    fname = joinpath(path, "Output.nc")
    Stats = K1D.NetCDFIO_Stats(fname, 1.0, parent(face_coord), parent(coord))

    # Solve the initial value problem for density profile
    ρ_profile = K1D.ρ_ivp(FT, kid_params, thermo_params)
    # Create the initial condition profiles
    init = map(coord -> K1D.init_1d_column(FT, common_params, kid_params, thermo_params, ρ_profile, coord.z), coord)

    # Create state vector and apply initial condition
    Y = CO.initialise_state(moisture, precip, init)

    # Create aux vector and apply initial condition
    aux = K1D.initialise_aux(
        FT,
        init,
        common_params,
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
    K1D.KiD_output(aux, 0.0)

    # Define callbacks for output
    callback_io = ODE.DiscreteCallback(K1D.condition_io, K1D.affect_io!; save_positions = (false, false))
    callbacks = ODE.CallbackSet(callback_io)

    # Collect all the tendencies into rhs function for ODE solver
    # based on model choices for the solved equations
    ode_rhs! = K1D.make_rhs_function(moisture, precip)

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
        plot_animation(z_centers, solver, aux, moisture, precip, K1D, output = plot_folder)
        plot_timeheight(string("experiments/KiD_driver/", output_folder, "/Output.nc"), output = plot_folder)
    end

    return solver

end
