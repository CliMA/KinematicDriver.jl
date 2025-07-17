"""
    A driver for running the kinematic 1D simulations.

    Try: julia --project=test/ test/experiments/KiD_driver.jl --help
    to see the list of command line arguments.
"""

import KinematicDriver

# Get the parameter values for the simulation
kid_driver_path = joinpath(pkgdir(KinematicDriver), "test", "experiments", "KiD_driver")
include(joinpath(kid_driver_path, "parse_commandline.jl"))

if !(@isdefined config)
    config = parse_commandline()
end

ft_choice = config["FLOAT_TYPE"]
@assert ft_choice ∈ ("Float64", "Float32") "Invalid float type: $ft_choice"

FT = ft_choice == "Float64" ? Float64 : Float32

include(joinpath(kid_driver_path, "run_KiD_simulation.jl"))

run_KiD_simulation(FT, config)
