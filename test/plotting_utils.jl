"""
Plotting utilities
"""

import NCDatasets as NC
import CloudMicrophysics.PrecipitationSusceptibility as CMPS

ENV["GKSwstype"] = "nul"
import ClimaCorePlots, Plots
Plots.GRBackend()
using CairoMakie, Printf

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
    profs = ds.group["profiles"]

    z_plt = Array(profs["zc"])
    t_plt = Array(profs["t"])
    nt = length(t_plt)
    ti = Observable(1)  # time index

    kg_to_g = 1e3
    m⁻³_to_cm⁻³ = 1e-6
    q_tot = Array(profs["q_tot"]) .* kg_to_g
    q_liq = Array(profs["q_liq"]) .* kg_to_g
    q_ice = Array(profs["q_ice"]) .* kg_to_g
    q_rai = Array(profs["q_rai"]) .* kg_to_g
    q_sno = Array(profs["q_sno"]) .* kg_to_g
    N_aer = Array(profs["N_aer"]) .* m⁻³_to_cm⁻³
    N_liq = Array(profs["N_liq"]) .* m⁻³_to_cm⁻³
    N_rai = Array(profs["N_rai"]) .* m⁻³_to_cm⁻³
    SN_liq_act = Array(profs["SN_liq_act"]) .* m⁻³_to_cm⁻³

    close(ds)

    lims(x) = iszero(x) ? (nothing, nothing) : ((0, 1.1 * maximum(x)), nothing)

    fig = Figure(size = (1800, 1500))
    ax_q_tot = Axis(fig[1, 1]; xlabel = "q_tot [g/kg]", limits = lims(q_tot), ylabel = "z [m]")
    lines!(ax_q_tot, @lift(q_tot[:, $ti]), z_plt)
    ax_q_liq = Axis(fig[1, 2]; xlabel = "q_liq [g/kg]", limits = lims(q_liq))
    lines!(ax_q_liq, @lift(q_liq[:, $ti]), z_plt)
    ax_N_liq = Axis(fig[1, 3]; xlabel = "N_liq [1/cm^3]", limits = lims(N_liq))
    lines!(ax_N_liq, @lift(N_liq[:, $ti]), z_plt)
    ax_q_rai = Axis(fig[2, 1]; xlabel = "q_rai [g/kg]", limits = lims(q_rai), ylabel = "z [m]")
    lines!(ax_q_rai, @lift(q_rai[:, $ti]), z_plt)
    ax_N_rai = Axis(fig[2, 2]; xlabel = "N_rai [1/cm^3]", limits = lims(N_rai))
    lines!(ax_N_rai, @lift(N_rai[:, $ti]), z_plt)
    ax_q_ice = Axis(fig[2, 3]; xlabel = "q_ice [g/kg]", limits = lims(q_ice))
    lines!(ax_q_ice, @lift(q_ice[:, $ti]), z_plt)
    ax_q_sno = Axis(fig[3, 1]; xlabel = "q_sno [g/kg]", limits = lims(q_sno), ylabel = "z [m]")
    lines!(ax_q_sno, @lift(q_sno[:, $ti]), z_plt)
    ax_SN_liq_act = Axis(fig[3, 2]; xlabel = "SN_liq_act [1/cm^3]", limits = lims(SN_liq_act))
    lines!(ax_SN_liq_act, @lift(SN_liq_act[:, $ti]), z_plt)
    axs = contents(fig[:, :])
    Label(fig[0, :], @lift(@sprintf("time = %d [s]", t_plt[$ti])), fontsize = 30)
    linkyaxes!(axs...)

    record(fig, joinpath(path, "animation.mp4"), 1:nt) do i
        ti[] = i
    end
end

function plot_profiles_in_time(nc_data_file; output = "output", n = 10)

    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    profs = ds.group["profiles"]

    z_plt = Array(profs["zc"])
    t_plt = Array(profs["t"])
    nt = length(t_plt)

    # Select n evenly distributed time indices
    Δi = round(Int, (nt - 2) / n)
    time_indices = [1; Δi:Δi:(nt - Δi); nt]  # explicitly include first and last indices
    ni = length(time_indices)
    selected_times = t_plt[time_indices]

    kg_to_g = 1e3
    m⁻³_to_cm⁻³ = 1e-6
    q_tot = profs["q_tot"][:, time_indices] .* kg_to_g
    q_liq = profs["q_liq"][:, time_indices] .* kg_to_g
    q_ice = profs["q_ice"][:, time_indices] .* kg_to_g
    q_rai = profs["q_rai"][:, time_indices] .* kg_to_g
    q_sno = profs["q_sno"][:, time_indices] .* kg_to_g
    N_liq = profs["N_liq"][:, time_indices] .* m⁻³_to_cm⁻³
    N_rai = profs["N_rai"][:, time_indices] .* m⁻³_to_cm⁻³
    SN_liq_act = profs["SN_liq_act"][:, time_indices] .* m⁻³_to_cm⁻³

    close(ds)

    # Create colormap for time progression
    # colormap = :viridis
    colormap = cgrad(:darkrainbow, (0:ni) / ni, categorical = true)
    # colors = get(colorschemes[colormap], range(0, 1, length = n))
    ymax = z_plt[end] + (z_plt[end] - z_plt[end - 1])  # extend y-axis a bit beyond the last point
    lims(x) = iszero(x) ? (nothing, (0, ymax)) : ((0, 1.1 * maximum(x)), (0, ymax))

    fig = Figure(size = (2000, 1500))

    # Create axes
    ax_q_tot = Axis(fig[1, 1]; xlabel = "q_tot [g/kg]", limits = lims(q_tot), ylabel = "z [m]")
    ax_q_liq = Axis(fig[1, 2]; xlabel = "q_liq [g/kg]", limits = lims(q_liq))
    ax_N_liq = Axis(fig[1, 3]; xlabel = "N_liq [1/cm³]", limits = lims(N_liq))
    ax_q_rai = Axis(fig[2, 1]; xlabel = "q_rai [g/kg]", limits = lims(q_rai), ylabel = "z [m]")
    ax_N_rai = Axis(fig[2, 2]; xlabel = "N_rai [1/cm³]", limits = lims(N_rai))
    ax_q_ice = Axis(fig[2, 3]; xlabel = "q_ice [g/kg]", limits = lims(q_ice))
    ax_q_sno = Axis(fig[3, 1]; xlabel = "q_sno [g/kg]", limits = lims(q_sno), ylabel = "z [m]")
    ax_SN_liq_act = Axis(fig[3, 2]; xlabel = "SN_liq_act [1/cm³]", limits = lims(SN_liq_act))

    # Link y-axes
    axs = [ax_q_tot, ax_q_liq, ax_N_liq, ax_q_rai, ax_N_rai, ax_q_ice, ax_q_sno]
    linkyaxes!(axs...)

    # Plot lines for each selected time
    args = (; colormap, colorrange = (0, 1))
    for (i, t) in enumerate(selected_times)
        color = (i - 0.5) / ni  # Normalize color for colormap
        lines!(ax_q_tot, q_tot[:, i], z_plt; args..., color)
        lines!(ax_q_liq, q_liq[:, i], z_plt; args..., color)
        lines!(ax_N_liq, N_liq[:, i], z_plt; args..., color)
        lines!(ax_q_rai, q_rai[:, i], z_plt; args..., color)
        lines!(ax_N_rai, N_rai[:, i], z_plt; args..., color)
        lines!(ax_q_ice, q_ice[:, i], z_plt; args..., color)
        lines!(ax_q_sno, q_sno[:, i], z_plt; args..., color)
        lines!(ax_SN_liq_act, SN_liq_act[:, i], z_plt; args..., color)
    end

    # Add colorbar
    tloc = ((1:ni) .- 0.5) / ni  # locations offset by 0.5 to center label on the color
    tlab = begin
        # Canonicalize time labels
        can_times = @. canonicalize(Second(selected_times))
        # shorten time labels
        map(can_times) do t
            s = string(t)
            s == "empty period" && return "0 s"
            s = replace(s, "seconds" => "s", "minutes" => "m", "hours" => "h", "days" => "d")
            s = replace(s, "second" => "s", "minute" => "m", "hour" => "h", "day" => "d")
            s = replace(s, "milli" => "m", "micro" => "μ", "nano" => "n")
        end
    end
    Colorbar(fig[:, 4]; colormap, colorrange = (0, 1), width = 20, ticks = (tloc, tlab))

    # Add title
    Label(fig[0, :], "Vertical Profiles at Multiple Time Steps", fontsize = 24)

    # Save figure
    save(joinpath(path, "profiles_multitime.png"), fig)

    return fig
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

function make_tz_hm(i, j, t, z, var; imax = 3, title = "", kw...)
    fig = current_figure()
    ax = Axis(fig[i, j]; 
        xlabel = i == imax ? "time [s]" : "",
        ylabel = j == 1 ? "z [m]" : "",
        yticklabelsvisible = j == 1,
        xticklabelsvisible = i == 3,
        title,
    )
    hm = CairoMakie.heatmap!(ax, t, z, var; colormap = :BuPu, kw...)
    Makie.Colorbar(fig[i, j], hm;
        halign = 1, valign = 0,
        height = Makie.Relative(0.3), width = 10, tellwidth = false, 
        ticks = WilkinsonTicks(3; k_min = 3, k_max = 6),
    )
end

function plot_timeheight(nc_data_file; output = "output", mixed_phase = true, pysdm = false)
    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    if pysdm
        t = Array(ds["time"])
        z = Array(ds["height"])
        q_liq = Array(ds["cloud water mixing ratio"])
        q_rai = Array(ds["rain water mixing ratio"])
        q_ice = Array(ds["rain water mixing ratio"]) * FT(0)
        q_sno = Array(ds["rain water mixing ratio"]) * FT(0)
        q_vap = Array(ds["water_vapour_mixing_ratio"]) * 1e3
        q_tot = q_vap + q_liq
        N_aer = Array(ds["na"])
        N_liq = Array(ds["nc"])
        N_rai = Array(ds["nr"])
        SN_liq_act = transpose(Array(ds["activating"]))
    else
        data = ds.group["profiles"]
        t = Array(data["t"])
        z = Array(data["zc"])
        q_tot = permutedims(Array(data["q_tot"]))
        q_liq = permutedims(Array(data["q_liq"]))
        q_ice = permutedims(Array(data["q_ice"]))
        q_rai = permutedims(Array(data["q_rai"]))
        q_sno = permutedims(Array(data["q_sno"]))
        N_aer = permutedims(Array(data["N_aer"]))
        N_liq = permutedims(Array(data["N_liq"]))
        N_rai = permutedims(Array(data["N_rai"]))
        SN_liq_act = permutedims(Array(data["SN_liq_act"]))
    end

    kg_to_g = 1e3
    m⁻³_to_cm⁻³ = 1e-6

    if mixed_phase
        fig = Figure(size = (1200, 1200))
        make_tz_hm(1, 1, t, z, q_tot .* kg_to_g; colorrange = (8, 15), title = "q_tot [g/kg]")
        make_tz_hm(1, 2, t, z, q_liq .* kg_to_g; colorrange = (0, 1), title = "q_liq [g/kg]")
        make_tz_hm(1, 3, t, z, q_ice .* kg_to_g; colorrange = (0, 0.25), title = "q_ice [g/kg]")
        make_tz_hm(2, 1, t, z, q_rai .* kg_to_g; colorrange = (0, 0.25), title = "q_rai [g/kg]")
        make_tz_hm(2, 2, t, z, q_sno .* kg_to_g; title = "q_sno [g/kg]")
        make_tz_hm(2, 3, t, z, N_aer .* m⁻³_to_cm⁻³; colorrange = (0, 100), title = "N_aer [1/cm³]")
        make_tz_hm(3, 1, t, z, N_liq .* m⁻³_to_cm⁻³; colorrange = (0, 50), title = "N_liq [1/cm³]")
        make_tz_hm(3, 2, t, z, N_rai .* m⁻³_to_cm⁻³; colorrange = (0, 1), title = "N_rai [1/cm³]")
        make_tz_hm(3, 3, t, z, SN_liq_act .* m⁻³_to_cm⁻³; title = "Activation [1/cm³/s]")

    else
        fig = Figure(size = (1200, 600))
        make_tz_hm(1, 1, t, z, q_tot .* kg_to_g; colorrange = (8, 15), title = "q_tot [g/kg]")
        make_tz_hm(1, 2, t, z, q_liq .* kg_to_g; colorrange = (0, 1), title = "q_liq [g/kg]")
        make_tz_hm(1, 3, t, z, q_rai .* kg_to_g; colorrange = (0, 0.25), title = "q_rai [g/kg]")
        make_tz_hm(2, 1, t, z, N_aer .* m⁻³_to_cm⁻³; colorrange = (0, 100), title = "N_aer [1/cm³]")
        make_tz_hm(2, 2, t, z, N_liq .* m⁻³_to_cm⁻³; colorrange = (0, 50), title = "N_liq [1/cm³]")
        make_tz_hm(2, 3, t, z, N_rai .* m⁻³_to_cm⁻³; colorrange = (0, 1), title = "N_rai [1/cm³]")
        make_tz_hm(3, 1, t, z, SN_liq_act .* m⁻³_to_cm⁻³; title = "Activation [1/cm³/s]")
    end
    fig

    axs = filter(x -> x isa Axis, contents(fig[:, :]))
    linkaxes!(axs...)
    save(joinpath(path, "timeheight.png"), fig)
    nothing
end

function plot_timeheight_no_ice_snow(nc_data_file; output = "output", pysdm = false)
    path = joinpath(@__DIR__, output)
    mkpath(path)

    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    if pysdm
        t = Array(ds["time"])
        z = Array(ds["height"])
        q_liq = Array(ds["cloud water mixing ratio"]) / 1e3
        q_rai = Array(ds["rain water mixing ratio"]) / 1e3
        q_vap = Array(ds["water_vapour_mixing_ratio"])
        q_tot = q_vap + q_liq
    else
        t = Array(ds.group["profiles"]["t"])
        z = Array(ds.group["profiles"]["zc"])
        q_tot = permutedims(Array(ds.group["profiles"]["q_tot"]))
        q_liq = permutedims(Array(ds.group["profiles"]["q_liq"]))
        q_rai = permutedims(Array(ds.group["profiles"]["q_rai"]))
    end

    kg_to_g = 1e3
    fig = Figure(size = (1200, 300))
    make_tz_hm(1, 1, t, z, q_tot .* kg_to_g; title = "q_tot [g/kg]", imax = 1, colorrange = (0, 1))
    make_tz_hm(1, 2, t, z, q_liq .* kg_to_g; title = "q_liq [g/kg]", imax = 1, colorrange = (0, 1))
    make_tz_hm(1, 3, t, z, q_rai .* kg_to_g; title = "q_rai [g/kg]", imax = 1, colorrange = (0, 0.25))
    axs = filter(x -> x isa Axis, contents(fig[:, :]))
    linkaxes!(axs...)
    save(joinpath(path, "timeheight.png"), fig)
    nothing
end
