"""
    A driver for running the kinematic 1D simulations.

    Try: julia --project=test/ test/experiments/KiD_driver.jl --help
    to see the list of command line arguments.
"""

include("../../src/Kinematic1D.jl")
include("../create_parameters.jl")

import NCDatasets
import OrdinaryDiffEq
import UnPack
import ArgParse

import ClimaCore
import Thermodynamics
import CloudMicrophysics
import CLIMAParameters

const CC = ClimaCore
const TD = Thermodynamics
const CP = CLIMAParameters
const NC = NCDatasets
const ODE = OrdinaryDiffEq
const KD = Kinematic1D
const AP = ArgParse

const FT = Float64

# Get the parameter values for the simulation
include("parse_commandline.jl")
opts = parse_commandline()

# Equations to solve for mositure and precipitation variables
moisture_choice = opts["moisture_choice"]
precipitation_choice = opts["precipitation_choice"]
if moisture_choice == "EquilibriumMoisture"
    moisture = KD.EquilibriumMoisture()
elseif moisture_choice == "NonEquilibriumMoisture"
    moisture = KD.NonEquilibriumMoisture()
else
    error("Invalid moisture choice: $moisture_choice")
end
if precipitation_choice == "NoPrecipitation"
    precip = KD.NoPrecipitation()
elseif precipitation_choice == "Precipitation0M"
    precip = KD.Precipitation0M()
elseif precipitation_choice == "Precipitation1M"
    precip = KD.Precipitation1M()
else
    error("Invalid precipitation choice: $precipitation_choice")
end

# Initialize the timestepping struct
TS = KD.TimeStepping(FT(opts["dt"]), FT(opts["dt_output"]), FT(opts["t_end"]))

# Create the coordinates
space, face_space = KD.make_function_space(FT, opts["z_min"], opts["z_max"], opts["n_elem"])
coord = CC.Fields.coordinate_field(space)
face_coord = CC.Fields.coordinate_field(face_space)

# Initialize the netcdf output Stats struct
output_folder = string("Output_", moisture_choice, "_", precipitation_choice)
path = joinpath(@__DIR__, output_folder)
mkpath(path)
fname = joinpath(path, "Output.nc")
Stats = KD.NetCDFIO_Stats(fname, 1.0, vec(face_coord), vec(coord))

# Instantiate CliMA Parameters and overwrite the defaults
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
params = create_parameter_set(
    path,
    toml_dict,
    FT,
    opts["w1"],
    opts["t1"],
    Int(opts["precip_sources"]),
    Int(opts["precip_sinks"]),
)

# Solve the initial value problem for density profile
ρ_profile = KD.ρ_ivp(FT, params)
# Create the initial condition profiles
init = map(coord -> KD.init_1d_column(FT, params, ρ_profile, coord.z), coord)

# Create state vector and apply initial condition
Y = KD.initialise_state(moisture, precip, init)

# Create aux vector and apply initial condition
aux = KD.initialise_aux(FT, init, params, TS, Stats, face_space, moisture)

# Output the initial condition
KD.KiD_output(aux, 0.0)

# Define callbacks for output
callback_io = ODE.DiscreteCallback(KD.condition_io, KD.affect_io!; save_positions = (false, false))
callbacks = ODE.CallbackSet(callback_io)

# Collect all the tendencies into rhs function for ODE solver
# based on model choices for the solved equations
ode_rhs! = KD.make_rhs_function(moisture, precip)

# Solve the ODE operator
problem = ODE.ODEProblem(ode_rhs!, Y, (opts["t_ini"], opts["t_end"]), aux)
solver = ODE.solve(
    problem,
    ODE.SSPRK33(),
    dt = TS.dt,
    saveat = TS.dt_io,
    callback = callbacks,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);

# Some basic plots
if opts["plotting_flag"] == true

    include("../plotting_utils.jl")

    plot_folder = string("experiments/", output_folder, "/figures/")

    z_centers = parent(CC.Fields.coordinate_field(space))
    plot_final_aux_profiles(z_centers, aux, output = plot_folder)
    plot_animation(z_centers, solver, aux, moisture, precip, KD, output = plot_folder)
    plot_timeheight(string("experiments/", output_folder, "/Output.nc"))
end
