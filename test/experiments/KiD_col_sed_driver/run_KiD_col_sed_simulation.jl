import OrdinaryDiffEq as ODE
import ClimaCore as CC
import CLIMAParameters as CP
import CloudMicrophysics.Parameters as CMP
import Kinematic1D
import Kinematic1D.K1DModel as KID
import Kinematic1D.K1DColSedModel as KCS

include(joinpath(pkgdir(Kinematic1D), "test", "create_parameters.jl"))
include("plotting_utils.jl")

function run_KiD_col_sed_simulation(::Type{FT}, opts) where {FT}

    # Equations to solve for precipitation variables
    precipitation_choice = opts["precipitation_choice"]
    rain_formation_choice = opts["rain_formation_choice"]
    sedimentation_choice = opts["sedimentation_choice"]

    # Decide the output flder name based on options
    output_folder = string("Output_", precipitation_choice)
    if precipitation_choice in ["Precipitation1M", "Precipitation2M"]
        output_folder = output_folder * "_" * rain_formation_choice
        if sedimentation_choice == "Chen2022"
            output_folder = output_folder * "_Chen2022"
        end
    end
    path = joinpath(@__DIR__, output_folder)
    mkpath(path)

    # Overwrite the defaults parameters based on options
    default_toml_dict = CP.create_toml_dict(FT, dict_type = "alias")
    toml_dict = override_toml_dict(
        path,
        default_toml_dict,
        Int(opts["precip_sources"]),
        Int(opts["precip_sinks"]),
        FT(opts["prescribed_Nd"]),
    )
    # Create Kinematic1D model parameters
    # (some of the CloudMicrophysics.jl parameters structs are created later based on model choices)
    kid_params = create_kid_parameters(toml_dict)

    precip = KCS.get_precipitation_type(
        FT,
        precipitation_choice,
        rain_formation_choice,
        sedimentation_choice,
        toml_dict,
    )

    # Initialize the timestepping struct
    TS = KID.TimeStepping(FT(opts["dt"]), FT(opts["dt_output"]), FT(opts["t_end"]))

    # Create the coordinates
    space, face_space = KID.make_function_space(FT, FT(opts["z_min"]), FT(opts["z_max"]), opts["n_elem"])
    coord = CC.Fields.coordinate_field(space)
    face_coord = CC.Fields.coordinate_field(face_space)

    # Initialize the netcdf output Stats struct
    fname = joinpath(path, "Output.nc")
    Stats = KID.NetCDFIO_Stats(fname, 1.0, parent(face_coord), parent(coord), 
        output_profiles = Dict(
        :ρ => "density",
        :q_liq => "q_liq",
        :q_rai => "q_rai",
        :N_liq => "N_liq",
        :N_rai => "N_rai",
    ),)

    # Create the initial condition profiles
    init = map(coord -> KCS.init_1d_column(FT, opts["qt"], opts["prescribed_Nd"], opts["rhod"], coord.z), coord)

    # Create state vector and apply initial condition
    Y = KCS.initialise_state(precip, init)

    # Create aux vector and apply initial condition
    aux = KCS.initialise_aux(
        FT,
        init,
        kid_params,
        TS,
        Stats,
        face_space,
    )

    # Output the initial condition
    KID.KiD_output(aux, 0.0)

    # Define callbacks for output
    callback_io = ODE.DiscreteCallback(KID.condition_io, KID.affect_io!; save_positions = (false, false))
    callbacks = ODE.CallbackSet(callback_io)

    # Collect all the tendencies into rhs function for ODE solver
    # based on model choices for the solved equations
    ode_rhs! = KCS.make_rhs_function(precip)

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
    plot_folder = joinpath(pkgdir(Kinematic1D), string("test/experiments/KiD_col_sed_driver/", output_folder, "/figures/"))
    plot_timeheight(string(output_folder, "/Output.nc"), output = plot_folder)

    return solver
end

opts = Dict(
    "qt" => 1e-3,
    "prescribed_Nd" => 1e8,
    "rhod" => 1.0,
    "precipitation_choice" => "Precipitation1M",
    "rain_formation_choice" => "CliMA_1M",
    "sedimentation_choice" => "CliMA_1M",
    "precip_sources" => true,
    "precip_sinks" => true,
    "z_min" => 0.0,
    "z_max" => 3000.0,
    "n_elem" => 60,
    "dt" => 1.0,
    "dt_output" => 30.0,
    "t_ini" => 0.0,
    "t_end" => 3600.0,
)
run_KiD_col_sed_simulation(Float64, opts);
