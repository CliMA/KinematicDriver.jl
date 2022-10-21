using Plots

import Kinematic1D
const KID = Kinematic1D

include("../config.jl")

config = get_config()
u_names = collect(keys(config["prior"]["parameters"]))
u_values = collect([v.mean for v in values(config["prior"]["parameters"])])
model_settings = config["model"]
model_settings["p0"] = 99000
model_settings["Nd"] = 100.0 * 1e6

n_heights = config["model"]["n_elem"]
dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_heights
heights = collect(range(config["model"]["z_min"] + dz / 2, config["model"]["z_max"] - dz / 2, n_heights))
times = collect(config["model"]["t_calib"])
n_times = length(times)
variables = ["rl", "rr"]
n_vars = length(variables)

dir = "/Users/sajjadazimi/Postdoc/Code/PySDM/01-rain_shaft/data/02-p990/varying_rhow_Nd/"
rhow = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

obs_matrix = Dict()
G_matrix = Dict()
for rw in rhow
    file_data_vector = zeros(n_vars * n_heights * n_times)
    count = 0
    foreach(readdir(dir)) do file
        filename = dir * file
        if occursin("rhow=" * string(rw) * "_pc=100_", filename)
            file_data_vector = file_data_vector .+ KID.get_single_obs_vector(filename, variables, heights, times)
            count = count + 1
        end
    end
    file_data_vector = file_data_vector ./ count
    obs_matrix[rw] = reshape(file_data_vector, n_heights * n_vars, n_times)

    model_settings["w1"] = rw
    ode_sol, aux = KID.run_KiD(u_values, u_names, model_settings)
    G = KID.ODEsolution2Gvector(ode_sol, aux, variables, ones(n_vars))
    G_matrix[rw] = reshape(G, n_heights * n_vars, n_times)
end

p = Array{Plots.Plot}(undef, n_vars)
lc = Dict(
    0.5 => "blue",
    1.0 => "red",
    1.5 => "yellow",
    2.0 => "cyan",
    2.5 => "purple",
    3.0 => "green",
    3.5 => "magenta",
    4.0 => "brown",
)
p[1] = plot()
for rw in rhow
    p[1] = plot!(
        times ./ 60,
        sum(G_matrix[rw][1:n_heights, :], dims = 1)[:] .* dz,
        ylabel = "height [km]",
        linecolor = lc[rw],
    )

    p[1] = plot!(
        times ./ 60,
        sum(obs_matrix[rw][1:n_heights, :], dims = 1)[:] .* dz,
        linestyle = :dash,
        linecolor = lc[rw],
    )
end

if n_vars == 2
    p[2] = plot()
    for rw in rhow
        p[2] = plot!(
            times ./ 60,
            sum(G_matrix[rw][(n_heights + 1):(2 * n_heights), :], dims = 1)[:] .* dz,
            xlabel = "time [min]",
            ylabel = "height [km]",
            linecolor = lc[rw],
        )

        p[2] = plot!(
            times ./ 60,
            sum(obs_matrix[rw][(n_heights + 1):(2 * n_heights), :], dims = 1)[:] .* dz,
            linestyle = :dash,
            linecolor = lc[rw],
        )
    end
end

fig = plot(p..., layout = (n_vars, 1), legend = false)
Plots.png(fig, "rwp-lwp_default.png")
plot!()
