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
    elseif precip isa CO.PrecipitationP3
        q_rai_end = parent(aux.microph_variables.q_rai)
        q_liqonice_end = parent(aux.microph_variables.q_liqonice)
        q_rim_end = parent(aux.microph_variables.q_rim)
        B_rim_end = parent(aux.microph_variables.B_rim)
        N_ice_end = parent(aux.microph_variables.N_ice)
        N_liq_end = parent(aux.microph_variables.N_liq)
        N_rai_end = parent(aux.microph_variables.N_rai)
    else
        q_rai_end = q_tot_end .* 0.0
        q_sno_end = q_tot_end .* 0.0
    end

    if precip isa CO.Precipitation2M
        N_aer_end = parent(aux.microph_variables.N_aer)
        N_liq_end = parent(aux.microph_variables.N_liq)
        N_rai_end = parent(aux.microph_variables.N_rai)
    else
        N_aer_end = q_tot_end .* 0.0
        N_liq_end = q_tot_end .* 0.0
        N_rai_end = q_tot_end .* 0.0
    end

    p1 = Plots.plot(q_tot_end .* 1e3, z_centers, xlabel = "q_tot [g/kg]", ylabel = "z [m]")
    p2 = Plots.plot(q_liq_end .* 1e3, z_centers, xlabel = "q_liq [g/kg]", ylabel = "z [m]")
    p3 = Plots.plot(q_ice_end .* 1e3, z_centers, xlabel = "q_ice [g/kg]", ylabel = "z [m]")
    p4 = Plots.plot(T_end, z_centers, xlabel = "T [K]", ylabel = "z [m]")
    p5 = Plots.plot(q_rai_end .* 1e3, z_centers, xlabel = "q_rai [g/kg]", ylabel = "z [m]")

    if precip isa CO.PrecipitationP3
        p6 = Plots.plot(N_ice_end .* 1e-6, z_centers, xlabel = "N_ice [1/cm3]", ylabel = "z [m]")
        p11 = Plots.plot(N_ice_end .* 1e-6, z_centers, xlabel = "N_ice [1/cm3]", ylabel = "z [m]")
        p12 = Plots.plot(q_liqonice_end .* 1e3, z_centers, xlabel = "q_liqonice [g/kg]", ylabel = "z [m]")
        p13 = Plots.plot(q_rim_end .* 1e3, z_centers, xlabel = "q_rim [g/kg]", ylabel = "z [m]")
        p14 = Plots.plot(B_rim_end, z_centers, xlabel = "B_rim [-]", ylabel = "z [m]")
        p = Plots.plot(
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p11,
            p12,
            p13,
            p14,
            size = (1800.0, 1200.0),
            bottom_margin = 40.0 * Plots.PlotMeasures.px,
            left_margin = 80.0 * Plots.PlotMeasures.px,
            layout = (2, 5),
        )
        Plots.png(p, joinpath(path, "final_aux_profiles.png"))
    else
        p6 = Plots.plot(q_sno_end .* 1e3, z_centers, xlabel = "q_sno [g/kg]", ylabel = "z [m]")
        p8 = Plots.plot(N_aer_end .* 1e-6, z_centers, xlabel = "N_aer [1/cm3]", ylabel = "z [m]")
        p9 = Plots.plot(N_liq_end .* 1e-6, z_centers, xlabel = "N_liq [1/cm3]", ylabel = "z [m]")
        p10 = Plots.plot(N_rai_end .* 1e-6, z_centers, xlabel = "N_rai [1/cm3]", ylabel = "z [m]")

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
            p8,
            p9,
            p10,
            size = (1800.0, 1200.0),
            bottom_margin = 40.0 * Plots.PlotMeasures.px,
            left_margin = 80.0 * Plots.PlotMeasures.px,
        )
        Plots.png(p, joinpath(path, "final_aux_profiles.png"))
    end
end

function plot_animation_p3(z_centers, solver, aux, moisture, precip, K1D, output = plot_folder)

    path = joinpath(@__DIR__, output)
    mkpath(path)

    ρ = parent(aux.thermo_variables.ρ)
    nz = length(z_centers)
    nt = length(solver.u)
    q_tot = zeros(nz, nt)
    q_liq = zeros(nz, nt)
    q_ice = zeros(nz, nt)
    q_rai = zeros(nz, nt)
    q_sno = zeros(nz, nt)
    N_rai = zeros(nz, nt)
    N_liq = zeros(nz, nt)
    N_ice = zeros(nz, nt)
    ρq_tot = zeros(nz, nt)
    ρq_liq = zeros(nz, nt)
    ρq_ice = zeros(nz, nt)
    ρq_liqonice = zeros(nz, nt)
    ρq_rim = zeros(nz, nt)
    ρq_rai = zeros(nz, nt)
    B_rim = zeros(nz, nt)

    for (i, u) in enumerate(solver.u)


        ρq_tot[:, i] = parent(u.ρq_tot) .* 1e3
        ρq_liq[:, i] = parent(u.ρq_liq) .* 1e3
        ρq_ice[:, i] = parent(u.ρq_ice) .* 1e3
        ρq_liqonice[:, i] = parent(u.ρq_liqonice) .* 1e3
        ρq_rim[:, i] = parent(u.ρq_rim) .* 1e3
        ρq_rai[:, i] = parent(u.ρq_rai) .* 1e3
        B_rim[:, i] = parent(u.B_rim) .* 1e3
        N_rai[:, i] = parent(u.N_rai) ./ 1e6
        N_liq[:, i] = parent(u.N_liq) ./ 1e6
        N_ice[:, i] = parent(u.N_ice) ./ 1e6

    t_plt = Array(ds.group["profiles"]["t"])
    z_plt = Array(ds.group["profiles"]["zc"])
    q_tot = Array(ds.group["profiles"]["q_tot"])
    q_liq = Array(ds.group["profiles"]["q_liq"])
    q_ice = Array(ds.group["profiles"]["q_ice"])
    q_rai = Array(ds.group["profiles"]["q_rai"])
    q_sno = Array(ds.group["profiles"]["q_sno"])
    N_aer = Array(ds.group["profiles"]["N_aer"])
    N_liq = Array(ds.group["profiles"]["N_liq"])
    N_rai = Array(ds.group["profiles"]["N_rai"])

    function plot_data(data, data_label, max_val, title = "")
        return Plots.plot(
            data,
            z_plt,
            xlabel = data_label,
            ylabel = "z [m]",
            legend = false,
            title = title,
            titlefontsize = 30,
            xlim = [0, 1.1 * max_val],
        )
    end

    anim = Plots.@animate for i in 1:length(t_plt)

        title = "time = " * string(floor(Int, solver.t[i])) * " [s]"

        p1 = plot_data(ρq_tot[:, i], "ρq_tot [g/m3]", maximum(ρq_tot))
        p2 = plot_data(ρq_liq[:, i], "ρq_liq [g/m3]", maximum(ρq_liq))
        p3 = plot_data(ρq_ice[:, i], "ρq_ice [g/m3]", maximum(ρq_ice), title)
        p4 = plot_data(ρq_liqonice[:, i], "ρq_liqonice [g/m3]", maximum(ρq_liqonice))
        p5 = plot_data(ρq_rim[:, i], "ρq_rim [g/m3]", maximum(ρq_rim))
        p6 = plot_data(ρq_rai[:, i], "ρq_rai [g/m3]", maximum(ρq_rai))
        p7 = plot_data(B_rim[:, i], "B_rim [-]", maximum(B_rim))
        p8 = plot_data(N_liq[:, i], "N_liq [1/cm^3]", maximum(N_liq))
        p9 = plot_data(N_ice[:, i], "N_ice [1/cm^3]", maximum(N_ice))
        p10 = plot_data(N_rai[:, i], "N_rai [1/cm^3]", maximum(N_rai))
        plot(
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
            size = (1800.0, 1500.0),
            bottom_margin = 30.0 * Plots.PlotMeasures.px,
            left_margin = 30.0 * Plots.PlotMeasures.px,
            top_margin = 30.0 * Plots.PlotMeasures.px,
            right_margin = 30.0 * Plots.PlotMeasures.px,
            layout = (5, 2),
        )
    end

    Plots.mp4(anim, joinpath(path, "animation.mp4"), fps = 10)
end

function plot_timeheight(nc_data_file; output = "output", mixed_phase = true, pysdm = false)

    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    if pysdm
        t_plt = Array(ds["time"])
        z_plt = Array(ds["height"])
        q_liq_plt = transpose(Array(ds["cloud water mixing ratio"]))
        q_rai_plt = transpose(Array(ds["rain water mixing ratio"]))
        q_ice_plt = transpose(Array(ds["rain water mixing ratio"])) * FT(0)
        q_sno_plt = transpose(Array(ds["rain water mixing ratio"])) * FT(0)
        q_vap = transpose(Array(ds["water_vapour_mixing_ratio"])) * 1e3
        q_tot_plt = q_vap + q_liq_plt
        N_aer_plt = transpose(Array(ds["na"]))
        N_liq_plt = transpose(Array(ds["nc"]))
        N_rai_plt = transpose(Array(ds["nr"]))
    else
        t_plt = Array(ds.group["profiles"]["t"])
        z_plt = Array(ds.group["profiles"]["zc"])
        q_tot_plt = Array(ds.group["profiles"]["q_tot"])
        q_liq_plt = Array(ds.group["profiles"]["q_liq"])
        q_ice_plt = Array(ds.group["profiles"]["q_ice"])
        q_rai_plt = Array(ds.group["profiles"]["q_rai"])
        q_sno_plt = Array(ds.group["profiles"]["q_sno"])
        N_aer_plt = Array(ds.group["profiles"]["N_aer"])
        N_liq_plt = Array(ds.group["profiles"]["N_liq"])
        N_rai_plt = Array(ds.group["profiles"]["N_rai"])
    end
    #! format: off
    p1 = Plots.heatmap(t_plt, z_plt, q_tot_plt .* 1e3, title = "q_tot [g/kg]", xlabel = "time [s]", ylabel = "z [m]", color = :BuPu, clims=(8, 15))
    p2 = Plots.heatmap(t_plt, z_plt, q_liq_plt .* 1e3, title = "q_liq [g/kg]", xlabel = "time [s]", ylabel = "z [m]", color = :BuPu, clims=(0, 1))
    p3 = Plots.heatmap(t_plt, z_plt, q_ice_plt .* 1e3, title = "q_ice [g/kg]", xlabel = "time [s]", ylabel = "z [m]", color = :BuPu, clims=(0, 0.25))
    p4 = Plots.heatmap(t_plt, z_plt, q_rai_plt .* 1e3, title = "q_rai [g/kg]", xlabel = "time [s]", ylabel = "z [m]", color = :BuPu, clims=(0, 0.25))
    p5 = Plots.heatmap(t_plt, z_plt, q_sno_plt .* 1e3, title = "q_sno [g/kg]", xlabel = "time [s]", ylabel = "z [m]", color = :BuPu)
    p6 = Plots.heatmap(t_plt, z_plt, N_aer_plt .* 1e-6, title = "N_aer [1/cm^3]", xlabel = "time [s]", ylabel = "z [m]", color = :BuPu, clims=(0, 100))
    p7 = Plots.heatmap(t_plt, z_plt, N_liq_plt .* 1e-6, title = "N_liq [1/cm^3]", xlabel = "time [s]", ylabel = "z [m]", color = :BuPu, clims=(0, 50))
    p8 = Plots.heatmap(t_plt, z_plt, N_rai_plt .* 1e-6, title = "N_rai [1/cm^3]", xlabel = "time [s]", ylabel = "z [m]", color = :BuPu, clims=(0, 1))
    #! format: on
    if mixed_phase
        p = Plots.plot(
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            size = (1200.0, 1200.0),
            bottom_margin = 30.0 * Plots.PlotMeasures.px,
            left_margin = 30.0 * Plots.PlotMeasures.px,
            layout = (3, 3),
        )
    else
        p = Plots.plot(
            p1,
            p2,
            p4,
            p6,
            p7,
            p8,
            size = (1200.0, 600.0),
            bottom_margin = 30.0 * Plots.PlotMeasures.px,
            left_margin = 30.0 * Plots.PlotMeasures.px,
            layout = (2, 3),
        )
    end
    Plots.png(p, joinpath(path, "timeheight.png"))
end

function plot_timeheight_p3(nc_data_file, precip; output = "output")
    path = joinpath(@__DIR__, output)
    mkpath(path)
    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    t_plt = Array(ds.group["profiles"]["t"])
    z_plt = Array(ds.group["profiles"]["zc"])
    q_tot_plt = Array(ds.group["profiles"]["q_tot"])
    q_liq_plt = Array(ds.group["profiles"]["q_liq"])
    ρq_ice_plt = Array(ds.group["profiles"]["ρq_ice"])
    ρq_rim_plt = Array(ds.group["profiles"]["ρq_rim"])
    ρq_liqonice_plt = Array(ds.group["profiles"]["ρq_liqonice"])
    q_rai_plt = Array(ds.group["profiles"]["q_rai"])
    q_vap_plt = Array(ds.group["profiles"]["q_vap"])
    N_aer_plt = Array(ds.group["profiles"]["N_aer"])
    N_liq_plt = Array(ds.group["profiles"]["N_liq"])
    N_rai_plt = Array(ds.group["profiles"]["N_rai"])
    N_ice_plt = Array(ds.group["profiles"]["N_ice"])
    B_rim_plt = Array(ds.group["profiles"]["B_rim"])
    #! format: off
    p1 = Plots.heatmap(t_plt, z_plt, q_tot_plt .* 1e3, title = "q_tot [g/kg]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p2 = Plots.heatmap(t_plt, z_plt, q_liq_plt .* 1e3, title = "q_liq [g/kg]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p3 = Plots.heatmap(t_plt, z_plt, ρq_ice_plt .* 1e3, title = "ρq_ice [g/m3]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p4 = Plots.heatmap(t_plt, z_plt, q_rai_plt .* 1e3, title = "q_rai [g/kg]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p5 = Plots.heatmap(t_plt, z_plt, q_vap_plt .* 1e3, title = "q_vap [g/kg]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p9 = Plots.heatmap(t_plt, z_plt, ρq_rim_plt .* 1e3, title = "ρq_rim [g/m3]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p10 = Plots.heatmap(t_plt, z_plt, ρq_liqonice_plt .* 1e3, title = "ρq_liqonice [g/m3]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p6 = Plots.heatmap(t_plt, z_plt, N_aer_plt .* 1e-6, title = "N_aer [1/cm3]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis, clims=(0, 100))
    p7 = Plots.heatmap(t_plt, z_plt, N_liq_plt .* 1e-6, title = "N_liq [1/cm3]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p8 = Plots.heatmap(t_plt, z_plt, N_rai_plt .* 1e-6, title = "N_rai [1/cm3]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p11 = Plots.heatmap(t_plt, z_plt, N_ice_plt .* 1e-6, title = "N_ice [1/cm3]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)
    p12 = Plots.heatmap(t_plt, z_plt, B_rim_plt, title = "B_rim [-]", xlabel = "time [s]", ylabel = "z [m]", color = :viridis)

    #! format: on
    p = Plots.plot(
        p1,
        p2,
        p3,
        p4,
        p5,
        p6,
        p7,
        p8,
        p9,
        p10,
        p11,
        p12,
        size = (2000.0, 1500.0),
        bottom_margin = 30.0 * Plots.PlotMeasures.px,
        left_margin = 30.0 * Plots.PlotMeasures.px,
        layout = (4, 3),
    )
    Plots.png(p, joinpath(path, "timeheight.png"))
end

function plot_timeheight_no_ice_snow(nc_data_file; output = "output", pysdm = false)
    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    if pysdm
        t_plt = Array(ds["time"])
        z_plt = Array(ds["height"])
        q_liq_plt = transpose(Array(ds["cloud water mixing ratio"])) / 1e3
        q_rai_plt = transpose(Array(ds["rain water mixing ratio"])) / 1e3
        q_vap = transpose(Array(ds["water_vapour_mixing_ratio"]))
        q_tot_plt = q_vap + q_liq_plt
    else
        t_plt = Array(ds.group["profiles"]["t"])
        z_plt = Array(ds.group["profiles"]["zc"])
        q_tot_plt = Array(ds.group["profiles"]["q_tot"])
        q_liq_plt = Array(ds.group["profiles"]["q_liq"])
        q_rai_plt = Array(ds.group["profiles"]["q_rai"])
    end

    p1 = Plots.heatmap(
        t_plt,
        z_plt,
        q_tot_plt .* 1e3,
        title = "q_tot [g/kg]",
        xlabel = "time [s]",
        ylabel = "z [m]",
        color = :BuPu,
        clims = (0, 1),
    )
    p2 = Plots.heatmap(
        t_plt,
        z_plt,
        q_liq_plt .* 1e3,
        title = "q_liq [g/kg]",
        xlabel = "time [s]",
        ylabel = "z [m]",
        color = :BuPu,
        clims = (0, 1),
    )
    p3 = Plots.heatmap(
        t_plt,
        z_plt,
        q_rai_plt .* 1e3,
        title = "q_rai [g/kg]",
        xlabel = "time [s]",
        ylabel = "z [m]",
        color = :BuPu,
        clims = (0, 0.25),
    )
    p = Plots.plot(
        p1,
        p2,
        p3,
        size = (1200.0, 300.0),
        bottom_margin = 30.0 * Plots.PlotMeasures.px,
        left_margin = 30.0 * Plots.PlotMeasures.px,
        layout = (1, 3),
    )
    Plots.png(p, joinpath(path, "timeheight.png"))
end
