"""
   Interface to ODE callbacks for handling output.
"""

function condition_io(u, t, integrator)
    UnPack.@unpack TS, io_info = integrator.p
    UnPack.@unpack Stats = io_info
    TS.dt_io += TS.dt
    io_flag = false
    if TS.dt_io > Stats.output_interval
        TS.dt_io = 0
        io_flag = true
    end
    return io_flag || t ≈ 0 || t ≈ TS.t_max
end

function affect_io!(integrator)
    KiD_output(integrator.p, integrator.t)

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end
