"""
Plotting utilities
"""

import NCDatasets as NC

ENV["GKSwstype"] = "nul"
using ClimaCorePlots, Plots
Plots.GRBackend()

function plot_timeheight(nc_data_file; output = "output")
    path = joinpath(@__DIR__, output)
    mkpath(path)

    p = Array{Plots.Plot}(undef, 2)
    ds = NC.NCDataset(joinpath(@__DIR__, nc_data_file))
    t_plt = collect(ds.group["profiles"]["t"])
    z_plt = collect(ds.group["profiles"]["zc"])
    q_liq_plt = collect(ds.group["profiles"]["q_liq"])
    q_rai_plt = collect(ds.group["profiles"]["q_rai"])
    p[1] = Plots.heatmap(t_plt, z_plt, q_liq_plt .* 1e3, title = "q_liq [g/kg]", xlabel = "time [s]", ylabel = "z [m]")
    p[2] = Plots.heatmap(t_plt, z_plt, q_rai_plt .* 1e3, title = "q_rai [g/kg]", xlabel = "time [s]", ylabel = "z [m]")

    p = Plots.plot(p..., size = (900.0, 350.0), bottom_margin = 5Plots.mm, left_margin = 4Plots.mm)
    Plots.png(p, joinpath(path, "timeheight.png"))
end
