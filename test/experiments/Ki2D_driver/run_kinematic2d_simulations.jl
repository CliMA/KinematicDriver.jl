import OrdinaryDiffEq as ODE
import ClimaCore as CC
import CLIMAParameters as CP
import CloudMicrophysics.Parameters as CMP
import Kinematic1D
import Kinematic1D.K1DModel as K1D
import Kinematic1D.K2DModel as K2D

ENV["GKSwstype"] = "nul"
using ClimaCorePlots, Plots
Plots.GRBackend()

include(joinpath(pkgdir(Kinematic1D), "test", "create_parameters.jl"))

function run_K2D_simulation(::Type{FT}, opts) where {FT}

    # Equations to solve for mositure and precipitation variables
    moisture_choice = opts["moisture_choice"]
    prognostics_choice = opts["prognostic_vars"]
    precipitation_choice = opts["precipitation_choice"]
    rain_formation_choice = opts["rain_formation_choice"]
    sedimentation_choice = opts["sedimentation_choice"]

    # Decide the output flder name based on options
    output_folder = string("Output_", moisture_choice, "_", prognostics_choice, "_", precipitation_choice)
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
        w1 = FT(opts["w1"]),
        t1 = FT(opts["t1"]),
        p0 = FT(opts["p0"]),
        precip_sources = Int(opts["precip_sources"]),
        precip_sinks = Int(opts["precip_sinks"]),
        qtot_flux_correction = 1,
        prescribed_Nd = FT(opts["prescribed_Nd"]),
        r_dry = FT(opts["r_dry"]),
        std_dry = FT(opts["std_dry"]),
        κ = FT(opts["kappa"]),
    )
    # Create Thermodynamics.jl and Kinematic1D model parameters
    # (some of the CloudMicrophysics.jl parameters structs are created later based on model choices)
    thermo_params = create_thermodynamics_parameters(toml_dict)
    common_params = create_common_parameters(toml_dict)
    kid_params = create_kid_parameters(toml_dict)
    air_params = CMP.AirProperties(FT, toml_dict)
    activation_params = CMP.AerosolActivationParameters(FT, toml_dict)

    moisture = CO.get_moisture_type(FT, moisture_choice, toml_dict)
    precip = CO.get_precipitation_type(
        FT,
        precipitation_choice,
        toml_dict,
        rain_formation_choice = rain_formation_choice,
        sedimentation_choice = sedimentation_choice,
    )

    # Initialize the timestepping struct
    TS = CO.TimeStepping(FT(opts["dt"]), FT(opts["dt_output"]), FT(opts["t_end"]))

    # set up function space
    hv_center_space, hv_face_space = K2D.make_function_space(
        FT,
        xlim = (0, opts["domain_width"]),
        zlim = (0, opts["domain_height"]),
        helem = opts["nx"],
        velem = opts["nz"],
        npoly = opts["npoly"],
    )
    coords = CC.Fields.coordinate_field(hv_center_space)
    face_coords = CC.Fields.coordinate_field(hv_face_space)

    # Initialize the netcdf output Stats struct
    fname = joinpath(path, "Output.nc")
    Stats = CO.NetCDFIO_Stats(
        fname,
        1.0,
        parent(face_coords.z)[:, 1, 1, 1],
        parent(coords.z)[:, 1, 1, 1],
        x = parent(coords.x)[1, 1, 1, :],
    )

    # Solve the initial value problem for density profile
    ρ_profile = CO.ρ_ivp(FT, kid_params, thermo_params)
    # Create the initial condition profiles
    init =
        map(coord -> CO.initial_condition_1d(FT, common_params, kid_params, thermo_params, ρ_profile, coord.z), coords)

    # Create state vector and apply initial condition
    Y = CO.initialise_state(moisture, precip, init)

    # Create aux vector and apply initial condition
    aux = K2D.initialise_aux(
        FT,
        init,
        common_params,
        kid_params,
        thermo_params,
        air_params,
        activation_params,
        opts["domain_width"],
        opts["domain_height"],
        TS,
        Stats,
        hv_center_space,
        hv_face_space,
        moisture,
    )

    # Output the initial condition
    CO.simulation_output(aux, 0.0)

    # Define callbacks for output
    callback_io = ODE.DiscreteCallback(CO.condition_io, CO.affect_io!; save_positions = (false, false))
    callbacks = ODE.CallbackSet(callback_io)

    # Collect all the tendencies into rhs function for ODE solver
    # based on model choices for the solved equations
    ode_rhs! = K2D.make_rhs_function(moisture, precip)

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

    # post-processing
    qtot_max = maximum([maximum(solver.u[i].ρq_tot) for i in 1:length(solver.u)])
    qliq_max = maximum([maximum(solver.u[i].ρq_liq) for i in 1:length(solver.u)])
    qrai_max = maximum([maximum(solver.u[i].ρq_rai) for i in 1:length(solver.u)])
    anim = Plots.@animate for u in solver.u
        Plots.plot(u.ρq_tot, clim = (0, qtot_max))
    end
    Plots.mp4(anim, joinpath(path, "rhoqtot.mp4"), fps = 20)
    anim = Plots.@animate for u in solver.u
        Plots.plot(u.ρq_liq .* 1000.0, clim = (0, qliq_max .* 1000.0))
    end
    Plots.mp4(anim, joinpath(path, "rhoqliq.mp4"), fps = 20)
    anim = Plots.@animate for u in solver.u
        Plots.plot(u.ρq_rai .* 1000.0, clim = (0, qrai_max .* 1000.0))
    end
    Plots.mp4(anim, joinpath(path, "rhoqrai.mp4"), fps = 20)

    return solver
end

opts = Dict(
    "moisture_choice" => "NonEquilibriumMoisture",
    "prognostic_vars" => "RhodTQ",
    "precipitation_choice" => "Precipitation1M",
    "rain_formation_choice" => "CliMA_1M",
    "sedimentation_choice" => "CliMA_1M",
    "w1" => 3.0,
    "t1" => 1800.0,
    "p0" => 100000.0,
    "precip_sources" => true,
    "precip_sinks" => true,
    "prescribed_Nd" => 1e8,
    "r_dry" => 0.04 * 1e-6,
    "std_dry" => 1.4,
    "kappa" => 0.9,
    "domain_width" => 3000.0,
    "domain_height" => 3000.0,
    "nx" => 32,
    "nz" => 32,
    "npoly" => 4,
    "dt" => 1.0,
    "dt_output" => 120.0,
    "t_ini" => 0.0,
    "t_end" => 3600.0,
)
run_K2D_simulation(Float64, opts);
