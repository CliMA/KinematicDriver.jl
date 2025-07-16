"""
Plotting utilities
"""

import NCDatasets as NC
import CloudMicrophysics.PrecipitationSusceptibility as CMPS

ENV["GKSwstype"] = "nul"
import ClimaCorePlots, Plots
Plots.GRBackend()
using CairoMakie

function plot_initial_profiles_comparison(KM; sdm_case = "dry")
    sdm_data = load_sdm_data(sdm_case)
    path = joinpath(@__DIR__, "initial_condition_tests/output_init_profiles")
    mkpath(path)

    fig = Figure(; size = (1200, 750))
    z = vec(KM.z_centers)

    ax = Axis(fig[1, 1]; xlabel = "q_vap [g/kg]", ylabel = "z [m]")
    lines!(vec(KM.q_vap), z, label = "KM")
    lines!(vec(sdm_data.qv_sdm), sdm_data.z_sdm, label = "SDM")

    ax = Axis(fig[1, 2]; xlabel = "ρ [kg/m3]")
    lines!(vec(KM.ρ), z, label = "KM")
    lines!(vec(sdm_data.rho_sdm), sdm_data.z_sdm, label = "SDM")

    ax = Axis(fig[1, 3]; xlabel = "θ_dry [K]")
    lines!(vec(KM.θ_dry), z, label = "KM")
    lines!(vec(sdm_data.thetad_sdm), sdm_data.z_sdm, label = "SDM")

    ax = Axis(fig[2, 1]; xlabel = "T [K]", ylabel = "z [m]")
    lines!(vec(KM.T), z, label = "KM")
    lines!(vec(sdm_data.T_sdm), sdm_data.z_sdm, label = "SDM")

    ax = Axis(fig[2, 2]; xlabel = "p [hPa]")
    lines!(vec(KM.p ./ 100), z, label = "KM")
    lines!(vec(sdm_data.P_sdm ./ 100), sdm_data.z_sdm, label = "SDM")

    ax = Axis(fig[2, 3]; xlabel = "q_liq [g/kg]")
    lines!(vec(KM.q_liq), z, label = "KM")
    lines!(vec(sdm_data.ql_sdm), sdm_data.z_sdm, label = "SDM")

    axs = contents(fig[:, :])
    axislegend.(axs)
    linkyaxes!(axs...)

    save(joinpath(path, "$(sdm_case)_init_profile.png"), fig)
end

function plot_final_aux_profiles(z_centers, aux, precip; output = "output")

    vars = aux.microph_variables
    thermo_vars = aux.thermo_variables

    path = joinpath(@__DIR__, output)
    mkpath(path)

    T_end = vec(thermo_vars.T)
    q_tot_end = vec(vars.q_tot)
    ρ_end = vec(thermo_vars.ρ)

    q_liq_end = vec(vars.q_liq)
    q_ice_end = vec(vars.q_ice)

    # Allocate variables
    q_rai_end = zero(q_tot_end)
    q_sno_end = zero(q_tot_end)
    N_aer_end = zero(q_tot_end)
    N_liq_end = zero(q_tot_end)
    N_rai_end = zero(q_tot_end)

    if precip isa CO.Precipitation1M
        q_rai_end .= vec(vars.q_rai)
        q_sno_end .= vec(vars.q_sno)
    elseif precip isa CO.Precipitation2M
        q_rai_end .= vec(vars.q_rai)
        N_aer_end .= vec(vars.N_aer)
        N_liq_end .= vec(vars.N_liq)
        N_rai_end .= vec(vars.N_rai)
        args = (precip.rain_formation, q_liq_end, q_rai_end, ρ_end, N_liq_end)
        precip_sus_aut = CMPS.precipitation_susceptibility_autoconversion.(args...)
        precip_sus_acc = CMPS.precipitation_susceptibility_accretion.(args...)
        d_ln_pp_d_ln_q_liq_aut = getfield.(precip_sus_aut, :d_ln_pp_d_ln_q_liq)
        d_ln_pp_d_ln_q_rai_aut = getfield.(precip_sus_aut, :d_ln_pp_d_ln_q_rai)
        d_ln_pp_d_ln_q_liq_acc = getfield.(precip_sus_acc, :d_ln_pp_d_ln_q_liq)
        d_ln_pp_d_ln_q_rai_acc = getfield.(precip_sus_acc, :d_ln_pp_d_ln_q_rai)
    elseif precip isa CO.PrecipitationP3
        q_rai_end .= vec(vars.q_rai)
        N_ice_end .= vec(vars.N_ice)
        N_liq_end .= vec(vars.N_liq)
        N_rai_end .= vec(vars.N_rai)
        # additional variables for P3
        q_liqonice_end = vec(vars.q_liqonice)
        q_rim_end = vec(vars.q_rim)
        B_rim_end = vec(vars.B_rim)
    end

    kg_to_g = 1e3
    m⁻³_to_cm⁻³ = 1e-6

    fig = Figure(size = (1800, 1200))
    ax = Axis(fig[1, 1]; xlabel = "q_tot [g/kg]", ylabel = "z [m]")
    lines!(q_tot_end .* kg_to_g, z_centers)

    ax = Axis(fig[1, 2]; xlabel = "q_liq [g/kg]")
    lines!(q_liq_end .* kg_to_g, z_centers)

    ax = Axis(fig[1, 3]; xlabel = "q_ice [g/kg]")
    lines!(q_ice_end .* kg_to_g, z_centers)

    ax = Axis(fig[1, 4]; xlabel = "T [K]")
    lines!(T_end, z_centers)

    ax = Axis(fig[2, 1]; xlabel = "q_rai [g/kg]", ylabel = "z [m]")
    lines!(q_rai_end .* kg_to_g, z_centers)

    if precip isa CO.PrecipitationP3
        ax = Axis(fig[2, 2]; xlabel = "N_ice [1/cm³]")
        lines!(N_ice_end .* m⁻³_to_cm⁻³, z_centers)

        ax = Axis(fig[2, 3]; xlabel = "q_liqonice [g/kg]")
        lines!(q_liqonice_end .* kg_to_g, z_centers)

        ax = Axis(fig[2, 4]; xlabel = "q_rim [g/kg]")
        lines!(q_rim_end .* kg_to_g, z_centers)

        ax = Axis(fig[3, 1]; xlabel = "B_rim [-]", ylabel = "z [m]")
        lines!(B_rim_end, z_centers)
    else
        ax = Axis(fig[2, 2]; xlabel = "q_sno [g/kg]")
        lines!(q_sno_end .* kg_to_g, z_centers)

        ax = Axis(fig[2, 3]; xlabel = "N_aer [1/cm³]")
        lines!(N_aer_end .* m⁻³_to_cm⁻³, z_centers)

        ax = Axis(fig[2, 4]; xlabel = "N_liq [1/cm³]")
        lines!(N_liq_end .* m⁻³_to_cm⁻³, z_centers)

        ax = Axis(fig[3, 1]; xlabel = "N_rai [1/cm³]", ylabel = "z [m]")
        lines!(N_rai_end .* m⁻³_to_cm⁻³, z_centers)

        ax = Axis(fig[3, 2]; xlabel = "precipitation susceptibility", ylabel = "z [m]")
        if precip isa CO.Precipitation2M
            lines!(ax, d_ln_pp_d_ln_q_liq_aut, z_centers, label = "aut, q_liq", color = :red)
            lines!(ax, d_ln_pp_d_ln_q_rai_aut, z_centers, label = "aut, q_rai", color = :brown)
            lines!(ax, d_ln_pp_d_ln_q_liq_acc, z_centers, label = "acc, q_liq", color = :blue)
            lines!(ax, d_ln_pp_d_ln_q_rai_acc, z_centers, label = "acc, q_rai", color = :green)
            axislegend(ax, position = :lb)
        end
    end
    axs = contents(fig[:, :])
    linkyaxes!(axs...)
    save(joinpath(path, "final_aux_profiles.png"), fig)
    nothing
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

    end

    function plot_data(data, data_label, max_val, title = "")
        return Plots.plot(
            data,
            z_centers,
            xlabel = data_label,
            ylabel = "z [m]",
            legend = false,
            title = title,
            titlefontsize = 30,
            xlim = [0, 1.1 * max_val],
        )
    end

    anim = Plots.@animate for i in 1:nt

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
        Plots.plot(
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

function plot_animation(nc_data_file; output = "output")

    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))

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

        title = "time = " * string(floor(Int, t_plt[i])) * " [s]"
        p1 = plot_data(q_tot[:, i] .* 1e3, "q_tot [g/kg]", maximum(q_tot))
        p2 = plot_data(q_liq[:, i] .* 1e3, "q_liq [g/kg]", maximum(q_liq), title)
        p3 = plot_data(N_liq[:, i] .* 1e6, "N_liq [1/cm^3]", maximum(N_liq))
        p4 = plot_data(q_rai[:, i] .* 1e3, "q_rai [g/kg]", maximum(q_rai))
        p5 = plot_data(N_rai[:, i] .* 1e6, "N_rai [1/cm^3]", maximum(N_rai))
        p6 = plot_data(q_ice[:, i] .* 1e3, "q_ice [g/kg]", maximum(q_ice))
        p7 = plot_data(q_sno[:, i] .* 1e3, "q_sno [g/kg]", maximum(q_sno))

        Plots.plot(
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            size = (1800.0, 1500.0),
            bottom_margin = 30.0 * Plots.PlotMeasures.px,
            left_margin = 30.0 * Plots.PlotMeasures.px,
            top_margin = 30.0 * Plots.PlotMeasures.px,
            right_margin = 30.0 * Plots.PlotMeasures.px,
            layout = (3, 3),
        )
    end

    Plots.mp4(anim, joinpath(path, "animation.mp4"), fps = 10)
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
