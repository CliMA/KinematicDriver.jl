using DelimitedFiles
using Interpolations
using Plots

data = readdlm("initial_1D_conditions.txt", skipstart=1)

input_T = reverse(data[:, 1])
#Td = reverse(data[:, 2])
input_qv = reverse(data[:, 3])
#Qsat = reverse(data[:, 4])
input_p = reverse(data[:, 5])
input_z = reverse(data[:, 10])

sounding_T = linear_interpolation(input_z, input_T; extrapolation_bc = Line());
sounding_qv = linear_interpolation(input_z, input_qv; extrapolation_bc = Line());
sounding_p = linear_interpolation(input_z, input_p; extrapolation_bc = Line());

#=
#Plotting
y = range(minimum(input_z, maximum(input_z), length=1000)
x = sounding_T.(y)
p =plot(x, y, label="Interpolated", linewidth=2)
scatter!(T, Altitude, label="Sounding", markersize=1, color=:red)
xlabel!("Temperature (K)")
ylabel!("Altitude")
savefig(p, "temperature_vs_altitude.png")
println("Saved plot to temperature_vs_altitude.png")
=#

# TO DO: make this look nicer and neater
#        make a branch to actually be doing this in a base of code
#        figure out how to use this to actually generate initial conditions
