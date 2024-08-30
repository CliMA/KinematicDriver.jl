"""
    A driver for running the kinematic 1D simulations.

    Try: julia --project=test/ test/experiments/KiD_driver.jl --help
    to see the list of command line arguments.
"""


# Get the parameter values for the simulation
include("parse_commandline.jl")

bound = (;
    ice_start = false,
    _magnitude = Float64(0.5),
    _q_flux = Float64(0.65e-3),
    _N_flux = Float64(1e5),
    _F_rim = Float64(0.8),
    _F_liq = Float64(0.5),
    _ρ_r_init = Float64(900),
)

# if !(@isdefined config)
config = parse_commandline()
config["precipitation_choice"] = "PrecipitationP3"
config["moisture_choice"] = "MoistureP3"
config["n_elem"] = 36
config["z_max"] = 3000
config["t_end"] = 75
config["w1"] = 0
config["dt"] = Float64(0.5)
config["p3_boundary_condition"] = bound
# end

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
