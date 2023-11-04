using Plots
import Kinematic1D.CalibrateCMP as KCP

include("run_box_simulation.jl")

model_settings = Dict(
    "precipitation_choice" => "Precipitation2M",
    "rain_formation_choice" => "SB2006",
    "precip_sources" => true,
    "precip_sinks" => true,
    "dt" => Float64(1),
    "dt_output" => Float64(30),
    "t_ini" => Float64(0),
    "t_end" => Float64(3600),
    "rho_air" => Float64(1.22),
)

FT = Float64
qt = 1.0
Nd = 50
k = 5
pysdm_file =
    "/Users/sajjadazimi/Postdoc/Results/04-PySDM_0D_BoxModel/data/qt=" *
    string(qt) *
    "_Nd=" *
    string(Nd) *
    "_k=" *
    string(k) *
    ".nc"
ref = KCP.read_pysdm_data(pysdm_file).variables

model_settings["init_q_liq"] = FT(ref["qc"][1])
model_settings["init_q_rai"] = FT(ref["qr"][1])
model_settings["init_N_liq"] = FT(ref["nc"][1])
model_settings["init_N_rai"] = FT(ref["nr"][1])

solver = run_box_simulation(FT, model_settings)
time = solver.t
data = vcat(solver.u'...)

plot(ref["t"] ./ 60, ref["qc"] .* 1000, label = "SDM", lw = 2, ls = :dash)
for p in parameterizations
    plot!(time ./ 60, data[:, 1] .* 1000, label = p, lw = 2)
end
p1 = plot!(xlabel = "time [m]", ylabel = "qc [g/kg]", left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)

plot(ref["t"] ./ 60, ref["qr"] .* 1000, label = "SDM", lw = 2, ls = :dash)
for p in parameterizations
    plot!(time ./ 60, data[:, 2] .* 1000, label = p, lw = 2)
end
p2 = plot!(xlabel = "time [m]", ylabel = "qr [g/kg]", left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)

plot(ref["t"] ./ 60, ref["nc"] ./ 10^6, label = "SDM", lw = 2, ls = :dash)
plot!(time ./ 60, data[:, 3] ./ 10^6, label = "SB2006", lw = 2)
p3 = plot!(xlabel = "time [m]", ylabel = "Nc [1/cm^3]", right_margin = 3Plots.mm)


plot(ref["t"] ./ 60, ref["nr"] ./ 10^6, label = "SDM", lw = 2, ls = :dash)
plot!(time ./ 60, data[:, 4] ./ 10^6, label = "SB2006", lw = 2)
p4 = plot!(xlabel = "time [m]", ylabel = "Nr [1/cm^3]", right_margin = 3Plots.mm)

plot(p1, p2, p3, p4, size = (800, 600))
