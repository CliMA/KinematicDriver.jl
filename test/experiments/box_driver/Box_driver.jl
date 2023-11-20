"""
    A driver for running box simulations with collisons only.

    Try: julia --project=test/ test/experiments/Box_driver.jl --help
    to see the list of command line arguments.
"""

using Plots

include("run_box_simulation.jl")

model_settings = Dict(
    "precipitation_choice" => "Precipitation2M",
    "rain_formation_choice" => "SB2006",
    "init_q_liq" => Float64(1e-3),
    "init_q_rai" => Float64(0),
    "init_N_liq" => Float64(1e8),
    "init_N_rai" => Float64(0),
    "precip_sources" => true,
    "precip_sinks" => true,
    "dt" => Float64(1),
    "dt_output" => Float64(30),
    "t_ini" => Float64(0),
    "t_end" => Float64(3600),
    "rho_air" => Float64(1.22),
    "microphys_params" => Dict(),
)

FT = Float64
solver = run_box_simulation(FT, model_settings)

time = solver.t
data = vcat(solver.u'...)

plot(time ./ 60, data[:, 1] .* 1000, label = "qc", lw = 2)
plot!(time ./ 60, data[:, 2] .* 1000, label = "qr", lw = 2)
p1 = plot!(
    xlabel = "time [m]",
    ylabel = "specific water content [g/kg]",
    legend = :right,
    left_margin = 3Plots.mm,
    bottom_margin = 3Plots.mm,
)

if model_settings["precipitation_choice"] == "Precipitation1M"
    plot(p1, size = (400, 300))
elseif model_settings["precipitation_choice"] == "Precipitation2M"
    plot(time ./ 60, data[:, 3] ./ 10^6, xlabel = "time [m]", ylabel = "Nc [1/cm^3]", label = "N_c", lw = 2)
    plot!(twinx(), time ./ 60, data[:, 4] ./ 10^6, ylabel = "Nr [1/cm^3]", label = "N_r", c = 2, lw = 2)
    p2 = plot!(legend = false, right_margin = 3Plots.mm)
    plot(p1, p2, size = (800, 300))
end
