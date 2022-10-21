using Plots

import Kinematic1D
const KID = Kinematic1D

include("../config.jl")

root_dir = "./output/CliMA_1M_9_params/"
files = [
    "parameters_Nd=10.jld2",
    "parameters_Nd=20.jld2",
    "parameters_Nd=50.jld2",
    "parameters_Nd=100.jld2",
    "parameters_Nd=200.jld2",
    "parameters_Nd=500.jld2",
    "parameters_Nd=1000.jld2",
]

Nd = [10, 20, 50, 100, 200, 500, 1000]

file_name = root_dir * files[1]
data = KID.load_data(file_name)
u_names = data.u_names

params = Dict()
for (i, file) in enumerate(files)
    file_name = root_dir * file
    data = KID.load_data(file_name)
    u_values = data.u

    for (j, name) in enumerate(u_names)
        if i == 1
            params[name] = []
        end
        params[name] = [params[name]; u_values[j]]
    end
end

p = Array{Plots.Plot}(undef, length(u_names))

for (i, name) in enumerate(u_names)
    p[i] = plot(Nd, params[name], xaxis = :log, label = false, xlabel = "Nd", ylabel = name)
end

fig =
    plot(p..., layout = (3, 3), legend = false, size = (1200, 800), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
Plots.png(fig, "params_vs_Nd.png")
plot!()
