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
