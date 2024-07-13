import OrdinaryDiffEq as ODE
using Plots

import ClimaCore as CC
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import KinematicDriver
import KinematicDriver.BoxModel as BX

include(joinpath(pkgdir(KinematicDriver), "test", "create_parameters.jl"))

function run_box_simulation(::Type{FT}, opts) where {FT}

    # Equations to solve for precipitation variables
    precipitation_choice = opts["precipitation_choice"]
    rain_formation_choice = opts["rain_formation_choice"]

    # Decide the output flder name based on options
    output_folder = string("Output_", precipitation_choice)
    if precipitation_choice in ["Precipitation1M", "Precipitation2M"]
        output_folder = output_folder * "_" * rain_formation_choice
    end
    path = joinpath(@__DIR__, output_folder)
    mkpath(path)

    # Overwrite the defaults parameters based on options
    default_toml_dict = CP.create_toml_dict(FT)
    toml_dict = override_toml_dict(
        path,
        default_toml_dict,
        precip_sources = 1,
        precip_sinks = 0,
        prescribed_Nd = FT(opts["Nd"]),
    )
    for (k, v) in opts["microphys_params"]
        toml_dict[k]["value"] = v
    end
    common_params = create_common_parameters(toml_dict)
    thermo_params = create_thermodynamics_parameters(toml_dict)
    air_params = CMP.AirProperties(toml_dict)
    activation_params = CMP.AerosolActivationParameters(toml_dict)

    moisture = CO.get_moisture_type("NonEquilibriumMoisture", toml_dict)
    precip = CO.get_precipitation_type(precipitation_choice, toml_dict; rain_formation_choice = rain_formation_choice)

    # Initialize the timestepping struct
    TS = CO.TimeStepping(FT(opts["dt"]), FT(opts["dt_output"]), FT(opts["t_end"]))

    # Create the 0D box domain (one z point)
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(0),
        CC.Geometry.ZPoint{FT}(1),
        boundary_names = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems = 1)
    space = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    coord = CC.Fields.coordinate_field(space)

    # Create the initial condition profiles
    init =
        map(coord -> CO.initial_condition_0d(FT, thermo_params, opts["qt"], opts["Nd"], opts["k"], opts["rhod"]), coord)

    # Create state vector and apply initial condition
    Y = CO.initialise_state(moisture, precip, init)

    # Create aux vector and apply initial condition
    aux = CO.initialise_aux(
        FT,
        init,
        common_params,
        thermo_params,
        air_params,
        activation_params,
        TS,
        0.0,
        moisture,
        precip,
    )

    # Collect all the tendencies into rhs function for ODE solver
    # based on model choices for the solved equations
    ode_rhs! = BX.make_rhs_function(moisture, precip)

    # Solve the ODE operator
    problem = ODE.ODEProblem(ode_rhs!, Y, (FT(opts["t_ini"]), FT(opts["t_end"])), aux)
    solver = ODE.solve(problem, ODE.SSPRK33(), dt = TS.dt, saveat = TS.dt_io)

    # plot results
    time = solver.t
    q_liq = []
    q_rai = []
    for u in solver.u
        q_liq = [q_liq; parent(u.ρq_liq) / (opts["rhod"] .+ parent(u.ρq_tot))]
        q_rai = [q_rai; parent(u.ρq_rai) / (opts["rhod"] .+ parent(u.ρq_tot))]
    end

    plot(time ./ 60, q_liq .* 1000, label = "qc", lw = 2)
    plot!(time ./ 60, q_rai .* 1000, label = "qr", lw = 2)
    p1 = plot!(
        xlabel = "time [m]",
        ylabel = "specific water content [g/kg]",
        legend = :right,
        left_margin = 3Plots.mm,
        bottom_margin = 3Plots.mm,
    )

    if opts["precipitation_choice"] == "Precipitation1M"
        plot(p1, size = (400, 300))
    elseif opts["precipitation_choice"] == "Precipitation2M"
        N_liq = []
        N_rai = []
        for u in solver.u
            N_liq = [N_liq; parent(u.N_liq)]
            N_rai = [N_rai; parent(u.N_rai)]
        end
        plot(time ./ 60, N_liq ./ 10^6, xlabel = "time [m]", ylabel = "Nc [1/cm^3]", label = "N_c", lw = 2)
        plot!(twinx(), time ./ 60, N_rai ./ 10^6, ylabel = "Nr [1/cm^3]", label = "N_r", c = 2, lw = 2)
        p2 = plot!(legend = false, right_margin = 3Plots.mm)
        plot(p1, p2, size = (800, 300))
    end

    savefig(path * "/result.pdf")

    return solver
end

opts = Dict(
    "qt" => 1e-3,
    "Nd" => 1e8,
    "k" => 2.0,
    "rhod" => 1.22,
    "precipitation_choice" => "Precipitation2M",
    "rain_formation_choice" => "SB2006",
    "dt" => 1.0,
    "dt_output" => 30.0,
    "t_ini" => 0.0,
    "t_end" => 3600.0,
    "microphys_params" => Dict(),
)
solver = run_box_simulation(Float64, opts);
