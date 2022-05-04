include("KiD.jl")

const FT = Float64

# Instantiate CliMA Parameters
struct AEPS <: CP.AbstractEarthParameterSet end
params = AEPS()

# Set up the computational domain and time step
z_min = FT(0)
z_max = FT(2e3)
n_elem = 256
Δt = 1.0
t_ini = 0.0
t_end = 10.0 * 60

# Updraft momentum flux terms and initial conditions
w1 = 2 # m/s * kg/m3
t1 = 600 # s
w_params = (w1 = w1, t1 = t1)

domain = Domains.IntervalDomain(
    Geometry.ZPoint{FT}(z_min),
    Geometry.ZPoint{FT}(z_max),
    boundary_names = (:bottom, :top),
)
mesh = Meshes.IntervalMesh(domain, nelems = n_elem)

space = Spaces.CenterFiniteDifferenceSpace(mesh)
face_space = Spaces.FaceFiniteDifferenceSpace(space)

coord = Fields.coordinate_field(space)
face_coord = Fields.coordinate_field(face_space)

# solve the initial value problem for density profile
ρ_profile = ρ_ivp(FT, params)
# create the initial condition profiles
init = map(coord -> init_1d_column(FT, params, ρ_profile, coord.z), coord)
w = Geometry.WVector.(zeros(FT, face_space))

# initialoze state and aux
# set initial condition
Y = Fields.FieldVector(; q_tot = init.q_tot)
aux = Fields.FieldVector(;
    ρ = init.ρ,
    θ_liq_ice = init.θ_liq_ice,
    T = init.T,
    q_liq = init.q_liq,
    q_ice = init.q_ice,
    w = w,
    w_params = w_params,
    params = params,
)

# Solve the ODE operator
problem = ODEProblem(rhs!, Y, (t_ini, t_end), aux)
solver = solve(
    problem,
    SSPRK33(),
    dt = Δt,
    saveat = 10 * Δt,
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

z_centers = parent(Fields.coordinate_field(space))

anim = Plots.@animate for u in solver.u
    q_tot = parent(u.q_tot)
    Plots.plot(q_tot, z_centers)
end
Plots.mp4(anim, joinpath(path, "KM_qt.mp4"), fps = 10)

θ_liq_ice_end = parent(aux.θ_liq_ice)
T_end = parent(aux.T)
q_liq_end = parent(aux.q_liq)
q_ice_end = parent(aux.q_ice)
q_tot_end = parent(solver.u[end].q_tot)
Plots.png(
    Plots.plot(θ_liq_ice_end, z_centers),
    joinpath(path, "KM_θ_end.png"),
)
Plots.png(
    Plots.plot(q_tot_end, z_centers),
    joinpath(path, "KM_qt_end.png"),
)
Plots.png(
    Plots.plot(q_liq_end, z_centers),
    joinpath(path, "KM_ql_end.png"),
)
Plots.png(
    Plots.plot(q_ice_end, z_centers),
    joinpath(path, "KM_qi_end.png"),
)
Plots.png(
    Plots.plot(T_end, z_centers),
    joinpath(path, "KM_T_end.png"),
)
