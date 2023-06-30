using ClimaCorePlots, Plots
Plots.GRBackend()

using FiniteDiff

import CloudMicrophysics.Microphysics2M as CM2

include("../experiments/KiD_driver/parse_commandline.jl")

include("../experiments/KiD_driver/run_KiD_simulation.jl")

path = joinpath(@__DIR__, "Output_test_precip_production")
mkpath(path)

const FT = Float64

"""
This is a very simple metric of precipitation production that is calculated
by summing the ρq_rai at a height corresponding to the cloud base over time.
This ignores effects due to varying terminal velocity, and does not have the
correct units, but is sufficient for testing susceptibility since it is
approximately proportional to the precipitation production.
"""
function total_precipitation_production(Nd)
    opts = parse_commandline()

    opts["plotting_flag"] = false
    opts["precipitation_choice"] = "Precipitation2M"
    opts["rain_formation_scheme_choice"] = "SB2006"
    opts["prescribed_Nd"] = Nd

    cloud_base = 75

    println(opts["prescribed_Nd"])

    solver = run_KiD_simulation(FT, opts)

    precip = sum([parent(u.ρq_rai)[cloud_base] for u in solver.u])

    println("Precipitation production: ", precip)

    return precip
end

Nd_range = FT.(range(1e5, 1e8, 10))

# p = Plots.plot(Nd_range, total_precipitation_production.(Nd_range), yscale = :log10,
#                 xlabel = "Prescribed Nd", ylabel = "Precipitation Production Metric")
# Plots.png(p, joinpath(path, "precip_production.png"))

precip_sus =
    Nd_range ./ total_precipitation_production.(Nd_range) .*
    FiniteDiff.finite_difference_derivative.(total_precipitation_production, Nd_range)

p = Plots.plot(
    Nd_range,
    precip_sus,
    xlabel = "Prescribed Nd",
    ylabel = "Precipitation Susceptibility",
    label = "Observed (KiD)",
)

expected_precip_sus = fill!(similar(Nd_range), -2)
Plots.plot!(Nd_range, expected_precip_sus, label = "Expected (Glassmeier & Lohmann)")

Plots.png(p, joinpath(path, "precip_susceptibility.png"))
