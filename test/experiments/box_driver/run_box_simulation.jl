import OrdinaryDiffEq as ODE

import CLIMAParameters as CP
import Kinematic1D.BoxModel as BX

function run_box_simulation(::Type{FT}, model_settings) where {FT}

    # Equations to solve for precipitation variables
    precipitation_choice = model_settings["precipitation_choice"]
    rain_formation_choice = model_settings["rain_formation_choice"]

    # Instantiate CliMA Parameters and overwrite the defaults
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    for (k, v) in model_settings["microphys_params"]
        toml_dict[k]["value"] = v
    end
    box_params = BX.Parameters.BoxModelParameters(
        FT(model_settings["rho_air"]),
        Int(model_settings["precip_sources"]),
        Int(model_settings["precip_sinks"]),
        FT(model_settings["init_N_liq"]),
    )

    precip = CO.get_precipitation_type(FT, precipitation_choice, toml_dict; rain_formation_choice = rain_formation_choice)

    # Initialize the timestepping struct
    TS = BX.TimeStepping(FT(model_settings["dt"]), FT(model_settings["dt_output"]), FT(model_settings["t_end"]))

    # Create state vector and apply initial condition
    Y = BX.initialise_state(precip, model_settings)

    # Create aux vector and apply initial condition
    aux = BX.initialise_aux(box_params, TS)

    # Collect all the tendencies into rhs function for ODE solver
    # based on model choices for the solved equations
    ode_rhs! = BX.make_rhs_function(precip)

    # Solve the ODE operator
    problem = ODE.ODEProblem(ode_rhs!, Y, (FT(model_settings["t_ini"]), FT(model_settings["t_end"])), aux)
    solver = ODE.solve(problem, ODE.SSPRK33(), dt = TS.dt, saveat = TS.dt_io)

    return solver
end
