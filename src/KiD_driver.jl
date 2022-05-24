"""
This should be turned into a test file
"""

include("KiD.jl")

const FT = Float64

# Instantiate CliMA Parameters
struct EarthParameterSet{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end
CP.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
nt = (; MSLP = 100000.0)
params = EarthParameterSet(nt)

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
q_surf = 0.016
ρw0 = 0.0

# initialize the timestepping struct
TS = TimeStepping(FT(Δt), FT(Δt_output), FT(t_end))

# create the coordinates,
space, face_space = make_function_space(z_min, z_max, n_elem)

coord = CC.Fields.coordinate_field(space)
face_coord = CC.Fields.coordinate_field(face_space)

# initialize the netcdf output Stats struct
Stats = NetCDFIO_Stats("Output.nc", 1.0, vec(face_coord), vec(coord))

# solve the initial value problem for density profile
ρ_profile = ρ_ivp(FT, params)
# create the initial condition profiles
init = map(coord -> init_1d_column(FT, params, ρ_profile, coord.z), coord)

# create state vector and apply initial condition
Y = initialise_state(EquilibriumMoisture(), NoPrecipitation(), init)

# create aux vector and apply initial condition
aux = initialise_aux(init, params, w_params, q_surf, ρw0, TS, Stats, face_space)

# output the initial condition
KiD_output(aux, 0.0)

# Define callbacks for output
callback_io = ODE.DiscreteCallback(condition_io, affect_io!; save_positions = (false, false))
callbacks = ODE.CallbackSet(callback_io)

# collect all the tendencies into rhs function for ODE solver
# based on model choices for the solved equations
ode_rhs! = make_rhs_function(EquilibriumMoisture(), NoPrecipitation())

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

# TODO - delete below once we have NetCDF output
ENV["GKSwstype"] = "nul"
using ClimaCorePlots, Plots
Plots.GRBackend()

dir = "advect"
path = joinpath(@__DIR__, "output", dir)
mkpath(path)

z_centers = parent(CC.Fields.coordinate_field(space))

anim = Plots.@animate for u in solver.u
    ρq_tot = parent(u.ρq_tot)
    ρ = parent(aux.ρ)
    Plots.plot(ρq_tot ./ ρ, z_centers)
end
Plots.mp4(anim, joinpath(path, "KM_qt.mp4"), fps = 10)

ρ = parent(aux.ρ)
θ_liq_ice_end = parent(aux.θ_liq_ice)
T_end = parent(aux.T)
q_liq_end = parent(aux.q_liq)
q_ice_end = parent(aux.q_ice)
ρq_tot_end = parent(solver.u[end].ρq_tot)

Plots.png(Plots.plot(θ_liq_ice_end, z_centers), joinpath(path, "KM_θ_end.png"))
Plots.png(Plots.plot(ρq_tot_end ./ ρ, z_centers), joinpath(path, "KM_qt_end.png"))
Plots.png(Plots.plot(q_liq_end, z_centers), joinpath(path, "KM_ql_end.png"))
Plots.png(Plots.plot(q_ice_end, z_centers), joinpath(path, "KM_qi_end.png"))
Plots.png(Plots.plot(T_end, z_centers), joinpath(path, "KM_T_end.png"))
Plots.png(Plots.plot(ρ, z_centers), joinpath(path, "KM_ρ.png"))
