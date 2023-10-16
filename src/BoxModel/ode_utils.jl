"""
  Main buliding blocks for the ODE solver
"""

"""
   Struct for storing the:
    - timestepping timestep `dt`,
    - output timestep `dt_io` and
    - simulation time `t_max`.
"""

mutable struct TimeStepping{FT <: Real}
    dt::FT
    dt_io::FT
    t_max::FT
end

"""
   Interface to ODE solver. Returns the function needed to compute the
   right hand side of the solved ODE. The rhs is assembled via dispatch
   based on the precipitation types.
"""
function make_rhs_function(ps::AbstractPrecipitationStyle)
    function rhs!(dY, Y, aux, t)
        precompute_aux_precip!(ps, dY, Y, aux, t)
    end
    return rhs!
end

"""
    Interface to ODE solver. Returns the named tuple with
    ODE solver state variables. The state is created via dispatching
    on different precipitation types
"""
function initialise_state(sp::AbstractPrecipitationStyle, model_settings)
    error("initailisation not implemented for a given $sp")
end
function initialise_state(::Precipitation1M, model_settings)
    return [model_settings["init_q_liq"], model_settings["init_q_rai"]]
end
function initialise_state(::Precipitation2M, model_settings)
    return [
        model_settings["init_q_liq"],
        model_settings["init_q_rai"],
        model_settings["init_N_liq"],
        model_settings["init_N_rai"],
    ]
end

"""
   Interface to ODE solver. It initializes the auxiliary state.
   The auxiliary state is passed to ODE solver via the `p` parameter 
   of the ODEProblem.
"""
function initialise_aux(box_params, TS)

    return (; box_params = box_params, TS = TS)
end
