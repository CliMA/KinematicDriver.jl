function condition_io(u, t, integrator)
    UnPack.@unpack TS, Stats = integrator.p
    TS.dt_io += TS.dt
    io_flag = false
    if TS.dt_io > Stats.output_interval
        TS.dt_io = 0
        io_flag = true
    end
    return io_flag || t ≈ 0 || t ≈ TS.t_max
end

function affect_io!(integrator)
    UnPack.@unpack Stats, ρ, T = integrator.p
    t = integrator.t

    open_files(Stats)

    # TODO: remove `vars` hack that avoids
    # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    # opening/closing files every step should be okay. #removeVarsHack
    # TurbulenceConvection.io(sim) # #removeVarsHack
    write_simulation_time(Stats, t) # #removeVarsHack

    #write_field(Stats, "density", ρ, "profiles")
    #write_field(Stats, "temperature", T, "profiles")

    close_files(Stats)

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end
