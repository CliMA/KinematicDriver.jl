"""
    A driver for running the kinematic 1D simulations.

    Try: julia --project=test/ test/experiments/KiD_driver.jl --help
    to see the list of command line arguments.
"""


# Get the parameter values for the simulation
include("parse_commandline.jl")

#if !(@isdefined config)
    config = parse_commandline()
    config["precipitation_choice"] = "PrecipitationP3"
    config["moisture_choice"] = "MoistureP3"
    config["n_elem"] = 15
    config["z_max"] = 3000
    config["t_end"] = 50
    config["w1"] = 0
    config["rv_0"] = 0
    config["rv_1"] = 0
    config["rv_2"] = 0
#end

ft_choice = config["FLOAT_TYPE"]

if ft_choice == "Float64"
    const FT = Float64
elseif ft_choice == "Float32"
    const FT = Float32
else
    error("Invalid float type: $ft_choice")
end

include("run_KiD_simulation.jl")

run_KiD_simulation(FT, config)
