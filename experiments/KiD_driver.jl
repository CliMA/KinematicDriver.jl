"""
This should be turned into a test file
"""

include("../src/Kinematic1D.jl")

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
const KiD = Kinematic1D

const FT = Float64

# Instantiate CliMA Parameters and overwrite the defaults
params = KiD.params_overwrite

# Chose the equations to solve for mositure variables
# (EquilibriumMoisture or NonEquilibriumMoisture).
moisture = KiD.NonEquilibriumMoisture()
# Chose the equations to solve for precipitation variables
# (NoPrecipitation, Precipitation0M or Precipitation1M).
precip = KiD.Precipitation0M()

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
TS = KiD.TimeStepping(FT(Δt), FT(Δt_output), FT(t_end))

# create the coordinates,
space, face_space = KiD.make_function_space(FT, z_min, z_max, n_elem)

coord = CC.Fields.coordinate_field(space)
face_coord = CC.Fields.coordinate_field(face_space)

# initialize the netcdf output Stats struct
fname = joinpath(path, "Output.nc")
nc_outputs = ("density", "temperature", "pressure")
ts_outputs = ("TODO")
Stats = KiD.NetCDFIO_Stats(fname, 1.0, vec(face_coord), vec(coord), nc_outputs)

#solve the initial value problem for density profile
ρ_profile = KiD.ρ_ivp(FT, params)
# create the initial condition profiles
init = map(coord -> KiD.init_1d_column(FT, params, ρ_profile, coord.z), coord)

# create state vector and apply initial condition
Y = KiD.initialise_state(moisture, precip, init)

# create aux vector and apply initial condition
aux = KiD.initialise_aux(FT, init, params, w_params, TS, Stats, nc_outputs, ts_outputs, face_space, moisture)

# # output the initial condition
KiD.KiD_output(aux, 0.0)

# Define callbacks for output
# callback_io = ODE.DiscreteCallback(KiD.condition_io, KiD.affect_io!; save_positions = (false, false))
# callbacks = ODE.CallbackSet(callback_io)

# # collect all the tendencies into rhs function for ODE solver
# # based on model choices for the solved equations
# ode_rhs! = KiD.make_rhs_function(moisture, precip)

# # Solve the ODE operator
# problem = ODE.ODEProblem(ode_rhs!, Y, (t_ini, t_end), aux)
# solver = ODE.solve(
#     problem,
#     ODE.SSPRK33(),
#     dt = Δt,
#     saveat = 10 * Δt,
#     callback = callbacks,
#     progress = true,
#     progress_message = (dt, u, p, t) -> t,
# );

# # TODO - delete below once we have NetCDF output
# ENV["GKSwstype"] = "nul"
# using ClimaCorePlots, Plots
# Plots.GRBackend()

# z_centers = parent(CC.Fields.coordinate_field(space))

# anim = Plots.@animate for u in solver.u
#     ρq_tot = parent(u.ρq_tot)
#     ρ = parent(aux.constants.ρ)
#     Plots.plot(ρq_tot ./ ρ, z_centers)
# end
# Plots.mp4(anim, joinpath(path, "KM_qt.mp4"), fps = 10)

# ρ = parent(aux.constants.ρ)
# θ_liq_ice_end = parent(aux.constants.θ_liq_ice)
# T_end = parent(aux.moisture_variables.T)
# q_liq_end = parent(aux.moisture_variables.q_liq)
# q_ice_end = parent(aux.moisture_variables.q_ice)
# ρq_tot_end = parent(solver.u[end].ρq_tot)

# Plots.png(Plots.plot(θ_liq_ice_end, z_centers), joinpath(path, "KM_θ_end.png"))
# Plots.png(Plots.plot(ρq_tot_end ./ ρ, z_centers), joinpath(path, "KM_qt_end.png"))
# Plots.png(Plots.plot(q_liq_end, z_centers), joinpath(path, "KM_ql_end.png"))
# Plots.png(Plots.plot(q_ice_end, z_centers), joinpath(path, "KM_qi_end.png"))
# Plots.png(Plots.plot(T_end, z_centers), joinpath(path, "KM_T_end.png"))
# Plots.png(Plots.plot(ρ, z_centers), joinpath(path, "KM_ρ.png"))