using CSV, DataFrames
using Plots
using Tables

q_rai_plt_chen = CSV.read("q_rai_plt_chen.csv", DataFrame, header= false)
q_rai_end_chen = CSV.read("q_rai_end_chen.csv", DataFrame,header= false)
q_rai_middle_chen = CSV.read("q_rai_middle_chen.csv", DataFrame, header=false)
t_plt_chen = CSV.read("t_plt_chen.csv", DataFrame,header=false)
z_plt_chen = CSV.read("z_plt_chen.csv", DataFrame,header=false)
z_centers_end_chen = CSV.read("z_centers_end_chen.csv", DataFrame,header=false)
z_centers_middle_chen = CSV.read("z_centers_middle_chen.csv", DataFrame,header=false)

q_rai_plt_original = CSV.read("q_rai_plt_original.csv", DataFrame,header= false)
q_rai_end_original = CSV.read("q_rai_end_original.csv", DataFrame,header= false)
q_rai_middle_original = CSV.read("q_rai_middle_original.csv", DataFrame,header=false)
t_plt_original = CSV.read("t_plt_original.csv", DataFrame,header=false)
z_plt_original = CSV.read("z_plt_original.csv", DataFrame,header=false)
z_centers_end_original = CSV.read("z_centers_end_original.csv", DataFrame,header=false)
z_centers_middle_original = CSV.read("z_centers_middle_original.csv", DataFrame,header=false)

Plots.plot(Matrix(q_rai_end_chen) .* 1e3, Matrix(z_centers_end_original), xlabel = "q_rai [g/kg]", ylabel = "z [m]", label = "Chen")
Plots.plot!(Matrix(q_rai_end_original) .* 1e3, Matrix(z_centers_end_original), xlabel = "q_rai [g/kg]", ylabel = "z [m]", label = "Ogura")
Plots.savefig("final_aux_profiles_comparison.png")

Plots.plot(Matrix(q_rai_middle_chen), Matrix(z_centers_middle_original), xlabel = "q_rai [g/kg]", ylabel = "z [m]", label = "Chen")
Plots.plot!(Matrix(q_rai_middle_original), Matrix(z_centers_middle_original), xlabel = "q_rai [g/kg]", ylabel = "z [m]", label = "Ogura")
Plots.savefig("middle_aux_profiles_comparison.png")

p1 = Plots.heatmap(vec(Matrix(t_plt_original)), vec(Matrix(z_plt_original)), Matrix(q_rai_plt_original) .* 1e3, title = "q_rai [g/kg] - Original", xlabel = "time [s]", ylabel = "z [m]", clims = (0, 0.08))
p2 = Plots.heatmap(vec(Matrix(t_plt_chen)),vec(Matrix(z_plt_chen)), Matrix(q_rai_plt_chen) .* 1e3, title = "q_rai [g/kg] - Chen", xlabel = "time [s]", ylabel = "z [m]", clims = (0, 0.08))
p3 = Plots.heatmap(vec(Matrix(t_plt_chen)), vec(Matrix(z_plt_chen)), (Matrix(q_rai_plt_chen) - Matrix(q_rai_plt_original)) .* 1e3, title = "q_rai [g/kg] - Difference (Chen - Original)", titlefontsize = 10, xlabel = "time [s]", ylabel = "z [m]", c = cgrad([:blue, :white, :red], [0.4, 0, -0.4]))
p = Plots.plot(
    p1,
    p2,
    p3,
    size = (1500.0, 500.0),
    bottom_margin = 30.0 * Plots.PlotMeasures.px,
    left_margin = 30.0 * Plots.PlotMeasures.px,
    layout = (1, 3),
)
Plots.png(p, "timeheight_comparison.png")