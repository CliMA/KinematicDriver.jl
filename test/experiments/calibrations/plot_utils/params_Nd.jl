using Plots

import Kinematic1D
const KID = Kinematic1D

include("../config.jl")

root_dir = "/Users/sajjadazimi/Postdoc/Results/02-Calibrations_01/09-CliMA1M_da0p5dv0p4/"
w_dirs = ["01-w2/", "02-w3/", "03-w4/"]
files = [
    "01-Nd10/output/parameters_optim.jld2",
    "02-Nd20/output/parameters_optim.jld2",
    "03-Nd50/output/parameters_optim.jld2",
    "04-Nd100/output/parameters_optim.jld2",
    "05-Nd200/output/parameters_optim.jld2",
    "06-Nd500/output/parameters_optim.jld2",
    "07-Nd1000/output/parameters_optim.jld2",
]

Nd = [10, 20, 50, 100, 200, 500, 1000]

file_name = root_dir * w_dirs[1] * files[1]
data = KID.load_data(file_name)
u_names = data.u_names

params = Dict()
for name in u_names
    params[name] = zeros(length(w_dirs), length(Nd))
end
for (k, wdir) in enumerate(w_dirs)
    for (i, file) in enumerate(files)
        file_name = root_dir * wdir * file
        data = KID.load_data(file_name).u_bests
        u_values = data.ϕ_optim

        for (j, name) in enumerate(u_names)
            params[name][k, i] = u_values[j]
        end
    end
end

p = Array{Plots.Plot}(undef, length(u_names))

for (i, name) in enumerate(u_names)
    p[i] = scatter()
    for j in 1:length(w_dirs)
        if name == "τ_acnv_rai"
            p[i] = scatter!(
                Nd,
                params[name][j, :],
                xaxis = :log,
                yaxis = :log,
                label = false,
                xlabel = "Nd",
                ylabel = name,
                ylims = [1e3, 1e8],
            )
        else
            p[i] = scatter!(
                Nd,
                params[name][j, :],
                xaxis = :log,
                label = false,
                xlabel = "Nd",
                ylabel = name,
                ylims = [0, 1.1 .* maximum(params[name])],
            )
        end
    end
end

fig = plot(p..., layout = (3, 3), legend = false, size = (1000, 600), left_margin = 5Plots.mm)
Plots.pdf(fig, "params_vs_Nd.pdf")
plot!()
