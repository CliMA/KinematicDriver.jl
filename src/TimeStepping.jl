"""
   Struct for storing the:
    - timestepping timestep `dt`,
    - output timestep `dt_io` and
    - simulation time `t_max`.
"""

mutable struct TimeStepping
    dt::Float64
    dt_io::Float64
    t_max::Float64
end
