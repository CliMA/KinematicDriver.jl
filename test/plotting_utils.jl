"""
Plotting utilities
"""

import NCDatasets as NC
import CloudMicrophysics.PrecipitationSusceptibility as CMPS

ENV["GKSwstype"] = "nul"
using ClimaCorePlots, Plots
Plots.GRBackend()

function plot_initial_profiles_comparison(KM; sdm_case = "dry")
    sdm_data = load_sdm_data(sdm_case)
    path = joinpath(@__DIR__, "initial_condition_tests/output_init_profiles")
    mkpath(path)

    fig_name = string(sdm_case, "_init_profile.png")

    p1 = Plots.plot(KM.q_vap, KM.z_centers, label = "KM", xlabel = "q_vap [g/kg]", ylabel = "z [m]")
    Plots.plot!(p1, sdm_data.qv_sdm, sdm_data.z_sdm, label = "SDM")

    p2 = Plots.plot(KM.ρ, KM.z_centers, label = "KM ρ_tot", xlabel = "ρ [kg/m3]", ylabel = "z [m]")
    Plots.plot!(p2, sdm_data.rho_sdm, sdm_data.z_sdm, label = "SDM ρ")

    p3 = Plots.plot(KM.θ_dry, KM.z_centers, label = "KM", xlabel = "θ_dry [K]", ylabel = "z [m]")
    Plots.plot!(p3, sdm_data.thetad_sdm, sdm_data.z_sdm, label = "SDM")

    p4 = Plots.plot(KM.T, KM.z_centers, label = "KM", xlabel = "T [K]", ylabel = "z [m]")
    Plots.plot!(p4, sdm_data.T_sdm, sdm_data.z_sdm, label = "SDM")

    p5 = Plots.plot(KM.p ./ 100, KM.z_centers, label = "KM", xlabel = "p [hPa]", ylabel = "z [m]")
    Plots.plot!(p5, sdm_data.P_sdm ./ 100, sdm_data.z_sdm, label = "SDM")

    p6 = Plots.plot(KM.q_liq, KM.z_centers, label = "KM", xlabel = "q_liq [g/kg]", ylabel = "z [m]")
    Plots.plot!(p6, sdm_data.ql_sdm, sdm_data.z_sdm, label = "SDM")

    p = Plots.plot(
        p1,
        p2,
        p3,
        p4,
        p5,
        p6,
        size = (1200.0, 750.0),
        bottom_margin = 40.0 * Plots.PlotMeasures.px,
        left_margin = 80.0 * Plots.PlotMeasures.px,
    )
    Plots.png(p, joinpath(path, fig_name))
end

function plot_final_aux_profiles(z_centers, aux, precip; output = "output")

    path = joinpath(@__DIR__, output)
    mkpath(path)

    T_end = parent(aux.thermo_variables.T)
    q_tot_end = parent(aux.microph_variables.q_tot)
    ρ_end = parent(aux.thermo_variables.ρ)

    q_liq_end = parent(aux.microph_variables.q_liq)
    q_ice_end = parent(aux.microph_variables.q_ice)

    if precip isa CO.Precipitation1M
        q_rai_end = parent(aux.microph_variables.q_rai)
        q_sno_end = parent(aux.microph_variables.q_sno)
    elseif precip isa CO.Precipitation2M
        q_rai_end = parent(aux.microph_variables.q_rai)
        q_sno_end = q_tot_end .* 0.0
    else
        q_rai_end = q_tot_end .* 0.0
        q_sno_end = q_tot_end .* 0.0
    end

    p1 = Plots.plot(q_tot_end .* 1e3, z_centers, xlabel = "q_tot [g/kg]", ylabel = "z [m]")
    p2 = Plots.plot(q_liq_end .* 1e3, z_centers, xlabel = "q_liq [g/kg]", ylabel = "z [m]")
    p3 = Plots.plot(q_ice_end .* 1e3, z_centers, xlabel = "q_ice [g/kg]", ylabel = "z [m]")
    p4 = Plots.plot(T_end, z_centers, xlabel = "T [K]", ylabel = "z [m]")
    p5 = Plots.plot(q_rai_end .* 1e3, z_centers, xlabel = "q_rai [g/kg]", ylabel = "z [m]")
    p6 = Plots.plot(q_sno_end .* 1e3, z_centers, xlabel = "q_sno [g/kg]", ylabel = "z [m]")

    p7 = Plots.plot(xlabel = "precipitation susceptibility", ylabel = "z [m]")
    if precip isa CO.Precipitation2M
        N_liq_end = parent(aux.microph_variables.N_liq)
        precip_sus_aut =
            CMPS.precipitation_susceptibility_autoconversion.(
                Ref(precip.rain_formation),
                q_liq_end,
                q_rai_end,
                ρ_end,
                N_liq_end,
            )
        precip_sus_acc =
            CMPS.precipitation_susceptibility_accretion.(
                Ref(precip.rain_formation),
                q_liq_end,
                q_rai_end,
                ρ_end,
                N_liq_end,
            )
        Plots.plot!([r.d_ln_pp_d_ln_q_liq for r in precip_sus_aut], z_centers, label = "aut, q_liq", color = :red)
        Plots.plot!([r.d_ln_pp_d_ln_q_rai for r in precip_sus_aut], z_centers, label = "aut, q_rai", color = :brown)
        Plots.plot!([r.d_ln_pp_d_ln_q_liq for r in precip_sus_acc], z_centers, label = "acc, q_liq", color = :blue)
        Plots.plot!([r.d_ln_pp_d_ln_q_rai for r in precip_sus_acc], z_centers, label = "acc, q_rai", color = :green)
        Plots.plot!(legend = :outerright)
    end

    p = Plots.plot(
        p1,
        p2,
        p3,
        p4,
        p5,
        p6,
        p7,
        size = (1200.0, 750.0),
        bottom_margin = 40.0 * Plots.PlotMeasures.px,
        left_margin = 80.0 * Plots.PlotMeasures.px,
    )
    Plots.png(p, joinpath(path, "final_aux_profiles.png"))
end

function plot_animation(z_centers, solver, aux, moisture, precip, KiD; output = "output")

    path = joinpath(@__DIR__, output)
    mkpath(path)

    anim = Plots.@animate for u in solver.u

        ρ = parent(aux.thermo_variables.ρ)
        q_tot = parent(u.ρq_tot) ./ ρ .* 1e3

        if moisture isa CO.NonEquilibriumMoisture
            q_liq = parent(u.ρq_liq) ./ ρ .* 1e3
            q_ice = parent(u.ρq_ice) ./ ρ .* 1e3
        else
            q_liq = q_tot .* 0.0
            q_ice = q_tot .* 0.0
        end

        if precip isa CO.Precipitation1M
            q_rai = parent(u.ρq_rai) ./ ρ .* 1e3
            q_sno = parent(u.ρq_sno) ./ ρ .* 1e3
        elseif precip isa CO.Precipitation2M
            q_rai = parent(u.ρq_rai) ./ ρ .* 1e3
            q_sno = q_tot .* 0.0
        else
            q_rai = q_tot .* 0.0
            q_sno = q_tot .* 0.0
        end

        p1 = Plots.plot(q_tot, z_centers, xlabel = "q_tot [g/kg]", ylabel = "z [m]")
        p2 = Plots.plot(q_liq, z_centers, xlabel = "q_liq [g/kg]", ylabel = "z [m]")
        p3 = Plots.plot(q_ice, z_centers, xlabel = "q_ice [g/kg]", ylabel = "z [m]")
        p4 = Plots.plot(q_rai, z_centers, xlabel = "q_rai [g/kg]", ylabel = "z [m]")
        p5 = Plots.plot(q_sno, z_centers, xlabel = "q_sno [g/kg]", ylabel = "z [m]")

        p = Plots.plot(
            p1,
            p2,
            p3,
            p4,
            p5,
            size = (1500.0, 1000.0),
            bottom_margin = 30.0 * Plots.PlotMeasures.px,
            left_margin = 30.0 * Plots.PlotMeasures.px,
            layout = (2, 3),
        )
    end

    Plots.mp4(anim, joinpath(path, "animation.mp4"), fps = 10)
end

function plot_timeheight(nc_data_file; output = "output")
    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    t_plt = Array(ds.group["profiles"]["t"])
    z_plt = Array(ds.group["profiles"]["zc"])
    q_tot_plt = Array(ds.group["profiles"]["q_tot"])
    q_liq_plt = Array(ds.group["profiles"]["q_liq"])
    q_ice_plt = Array(ds.group["profiles"]["q_ice"])
    q_rai_plt = Array(ds.group["profiles"]["q_rai"])
    q_sno_plt = Array(ds.group["profiles"]["q_sno"])
    p1 = Plots.heatmap(t_plt, z_plt, q_tot_plt .* 1e3, title = "q_tot [g/kg]", xlabel = "time [s]", ylabel = "z [m]")
    p2 = Plots.heatmap(t_plt, z_plt, q_liq_plt .* 1e3, title = "q_liq [g/kg]", xlabel = "time [s]", ylabel = "z [m]")
    p3 = Plots.heatmap(t_plt, z_plt, q_ice_plt .* 1e3, title = "q_ice [g/kg]", xlabel = "time [s]", ylabel = "z [m]")
    p4 = Plots.heatmap(t_plt, z_plt, q_rai_plt .* 1e3, title = "q_rai [g/kg]", xlabel = "time [s]", ylabel = "z [m]")
    p5 = Plots.heatmap(t_plt, z_plt, q_sno_plt .* 1e3, title = "q_sno [g/kg]", xlabel = "time [s]", ylabel = "z [m]")
    p = Plots.plot(
        p1,
        p2,
        p3,
        p4,
        p5,
        size = (1500.0, 1000.0),
        bottom_margin = 30.0 * Plots.PlotMeasures.px,
        left_margin = 30.0 * Plots.PlotMeasures.px,
        layout = (2, 3),
    )
    Plots.png(p, joinpath(path, "timeheight.png"))
end
function plot_timeheight_no_ice_snow(nc_data_file; output = "output")
    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    t_plt = Array(ds.group["profiles"]["t"])
    z_plt = Array(ds.group["profiles"]["zc"])
    q_tot_plt = Array(ds.group["profiles"]["q_tot"])
    q_liq_plt = Array(ds.group["profiles"]["q_liq"])
    q_rai_plt = Array(ds.group["profiles"]["q_rai"])
    p1 = Plots.heatmap(t_plt, z_plt, q_tot_plt .* 1e3, title = "q_tot [g/kg]", xlabel = "time [s]", ylabel = "z [m]")
    p2 = Plots.heatmap(t_plt, z_plt, q_liq_plt .* 1e3, title = "q_liq [g/kg]", xlabel = "time [s]", ylabel = "z [m]")
    p3 = Plots.heatmap(t_plt, z_plt, q_rai_plt .* 1e3, title = "q_rai [g/kg]", xlabel = "time [s]", ylabel = "z [m]")
    p = Plots.plot(
        p1,
        p2,
        p3,
        size = (1000.0, 500.0),
        bottom_margin = 30.0 * Plots.PlotMeasures.px,
        left_margin = 30.0 * Plots.PlotMeasures.px,
        layout = (2, 3),
    )
    Plots.png(p, joinpath(path, "timeheight.png"))
end
