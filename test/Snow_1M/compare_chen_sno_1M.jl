using CSV, DataFrames
using Plots
using Tables

q_sno_plt_chen = CSV.read("q_sno_plt_chen.csv", DataFrame, header= false)
q_sno_end_chen = CSV.read("q_sno_end_chen.csv", DataFrame,header= false)
q_sno_middle_chen = CSV.read("q_sno_middle_chen.csv", DataFrame, header=false)
t_plt_sno_chen = CSV.read("t_plt_sno_chen.csv", DataFrame,header=false)
z_plt_sno_chen = CSV.read("z_plt_sno_chen.csv", DataFrame,header=false)
z_centers_sno_end_chen = CSV.read("z_centers_sno_end_chen.csv", DataFrame,header=false)
z_centers_sno_middle_chen = CSV.read("z_centers_sno_middle_chen.csv", DataFrame,header=false)

q_sno_plt_original = CSV.read("q_sno_plt_original.csv", DataFrame,header= false)
q_sno_end_original = CSV.read("q_sno_end_original.csv", DataFrame,header= false)
q_sno_middle_original = CSV.read("q_sno_middle_original.csv", DataFrame,header=false)
t_plt_sno_original = CSV.read("t_plt_sno_original.csv", DataFrame,header=false)
z_plt_sno_original = CSV.read("z_plt_sno_original.csv", DataFrame,header=false)
z_centers_sno_end_original = CSV.read("z_centers_sno_end_original.csv", DataFrame,header=false)
z_centers_sno_middle_original = CSV.read("z_centers_sno_middle_original.csv", DataFrame,header=false)

Plots.plot(Matrix(q_sno_end_chen) .* 1e3, Matrix(z_centers_sno_end_chen), xlabel = "q_sno [g/kg]", ylabel = "z [m]", label = "Chen")
Plots.plot!(Matrix(q_sno_end_original) .* 1e3, Matrix(z_centers_sno_end_original), xlabel = "q_sno [g/kg]", ylabel = "z [m]", label = "Ogura")
Plots.savefig("final_aux_profiles_sno_comparison.png")

Plots.plot(Matrix(q_sno_middle_chen), Matrix(z_centers_sno_middle_chen), xlabel = "q_sno [g/kg]", ylabel = "z [m]", label = "Chen")
Plots.plot!(Matrix(q_sno_middle_original), Matrix(z_centers_sno_middle_original), xlabel = "q_sno [g/kg]", ylabel = "z [m]", label = "Ogura")
Plots.savefig("middle_aux_profiles_sno_comparison.png")

p1 = Plots.heatmap(vec(Matrix(t_plt_sno_original)), vec(Matrix(z_plt_sno_original)), Matrix(q_sno_plt_original) .* 1e3, title = "q_sno [g/kg] - Original", xlabel = "time [s]", ylabel = "z [m]", clims = (0, 8))
p2 = Plots.heatmap(vec(Matrix(t_plt_sno_chen)),vec(Matrix(z_plt_sno_chen)), Matrix(q_sno_plt_chen) .* 1e3, title = "q_sno [g/kg] - Chen", xlabel = "time [s]", ylabel = "z [m]", clims = (0, 8))
p3 = Plots.heatmap(vec(Matrix(t_plt_sno_chen)), vec(Matrix(z_plt_sno_chen)), (Matrix(q_sno_plt_chen) - Matrix(q_sno_plt_original)) .* 1e3, title = "q_sno [g/kg] - Difference (Chen - Original)", titlefontsize = 10, xlabel = "time [s]", ylabel = "z [m]", c = cgrad([:blue, :white, :red], [0.05, 0, -0.05]))
p = Plots.plot(
    p1,
    p2,
    p3,
    size = (1500.0, 500.0),
    bottom_margin = 30.0 * Plots.PlotMeasures.px,
    left_margin = 30.0 * Plots.PlotMeasures.px,
    layout = (1, 3),
)
Plots.png(p, "timeheight_sno_comparison.png")