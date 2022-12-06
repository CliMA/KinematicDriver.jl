using Plots

import Kinematic1D
const KID = Kinematic1D

rootdir = "/Users/sajjadazimi/Postdoc/Code/Kinematic1D.jl/output/"
key = "parameters_EKP.jld2"
data = KID.load(rootdir * key)["u_bests"].cov_optim
u_names = [name[1:5] for name in KID.load(rootdir * key)["u_names"]]
num = size(data)[1]
cor_map = zeros(num, num)
for i in 1:num
    for j in num:-1:1
        if j > num - i
            cor_map[i, j] = NaN
        else
            cor_map[i, j] = data[i, num - j + 1]
        end
    end
end

fig = heatmap(
    cor_map,
    xticks = (collect(1:num), collect(u_names)[end:-1:2]),
    yticks = (collect(1:num), collect(u_names)[1:(end - 1)]),
    size = (800, 600),
)
Plots.png(fig, "correlation_full_scheme.png")
plot!()
