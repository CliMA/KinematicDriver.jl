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
    UnPack.@unpack Stats, ρ, T, θ_d, P, q_liq = integrator.p
    t = integrator.t

    open_files(Stats)

    # TODO: remove `vars` hack that avoids
    # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    # opening/closing files every step should be okay. #removeVarsHack
    # TurbulenceConvection.io(sim) # #removeVarsHack
    write_simulation_time(Stats, t) # #removeVarsHack

    write_field(Stats, "density", vec(ρ), "profiles")
    write_field(Stats, "temperature", vec(T), "profiles")
    write_field(Stats, "theta_dry", vec(θ_d), "profiles")
    write_field(Stats, "pressure", vec(P), "profiles")
    write_field(Stats, "q_liq", vec(q_liq), "profiles")
    write_field(Stats, "theta_ql", vec(θ_liq_ice), "profiles")

    close_files(Stats)

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end
