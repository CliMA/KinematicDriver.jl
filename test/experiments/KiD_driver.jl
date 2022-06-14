"""
This should be turned into a test file
"""

include("../../src/Kinematic1D.jl")

import NCDatasets
import OrdinaryDiffEq
import UnPack

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

const FT = Float64

# Instantiate CliMA Parameters and overwrite the defaults
params = KD.params_overwrite

# Chose the equations to solve for mositure variables
# (EquilibriumMoisture or NonEquilibriumMoisture).
moisture = KD.NonEquilibriumMoisture()
# Chose the equations to solve for precipitation variables
# (NoPrecipitation, Precipitation0M or Precipitation1M).
precip = KD.Precipitation0M()

# Output folder name TODO - make automatic
path = joinpath(@__DIR__, "Output_NonEquilibrium_Precipitation1M")
mkpath(path)

# Set up the computational domain and time step
z_min = FT(0)
z_max = FT(2e3)
n_elem = 256
Δt = 1.0
Δt_output = 10 * Δt
t_ini = 0.0
t_end = 10.0 * 60

# Updraft momentum flux terms, surface terms and initial conditions
w1 = 2 # m/s * kg/m3
t1 = 600 # s
w_params = (w1 = w1, t1 = t1)

# initialize the timestepping struct
TS = KD.TimeStepping(FT(Δt), FT(Δt_output), FT(t_end))

# create the coordinates,
space, face_space = KD.make_function_space(FT, z_min, z_max, n_elem)

coord = CC.Fields.coordinate_field(space)
face_coord = CC.Fields.coordinate_field(face_space)

# initialize the netcdf output Stats struct
fname = joinpath(path, "Output.nc")
Stats = KD.NetCDFIO_Stats(fname, 1.0, vec(face_coord), vec(coord))

# solve the initial value problem for density profile
ρ_profile = KD.ρ_ivp(FT, params)
# create the initial condition profiles
init = map(coord -> KD.init_1d_column(FT, params, ρ_profile, coord.z), coord)

# create state vector and apply initial condition
Y = KD.initialise_state(moisture, precip, init)

# create aux vector and apply initial condition
aux = KD.initialise_aux(FT, init, params, w_params, TS, Stats, face_space, moisture)

# output the initial condition
KD.KiD_output(aux, 0.0)

# Define callbacks for output
callback_io = ODE.DiscreteCallback(KD.condition_io, KD.affect_io!; save_positions = (false, false))
callbacks = ODE.CallbackSet(callback_io)

# collect all the tendencies into rhs function for ODE solver
# based on model choices for the solved equations
ode_rhs! = KD.make_rhs_function(moisture, precip)

# Solve the ODE operator
problem = ODE.ODEProblem(ode_rhs!, Y, (t_ini, t_end), aux)
solver = ODE.solve(
    problem,
    ODE.SSPRK33(),
    dt = Δt,
    saveat = 10 * Δt,
    callback = callbacks,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);

plotting_flag = true
if plotting_flag == true

    include("../plotting_utils.jl")

    z_centers = parent(CC.Fields.coordinate_field(space))
    plot_final_aux_profiles(z_centers, aux, output = "experiments/output/")
    plot_animation(z_centers, solver, aux, moisture, precip, KD, output = "experiments/output/")
end
