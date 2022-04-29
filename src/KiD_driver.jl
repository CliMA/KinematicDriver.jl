include("KiD.jl")

const FT = Float64

# Instantiate CliMA Parameters
struct AEPS <: APS end
params = AEPS()

# Set up the computational domain and time step
z_min = FT(0)
z_max = FT(2e3)
n_elem = 128
Δt = 10.0
t_ini = 0.0
t_end = 10.0

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

Yc = map(coord -> init_1d_column(FT, params, coord.z), coord)
w = Geometry.WVector.(zeros(FT, face_space)) #TODO - should be changing in time
Y = Fields.FieldVector(Yc = Yc, w = w)

# Solve the ODE operator
ODE_sys = (dY, Y, _, t) -> advection_tendency!(dY, Y, _, t, w_params)
problem = ODEProblem(ODE_sys, Y, (t_ini, t_end))
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

# anim = Plots.@animate for u in solver.u
#     θ = parent(u.Yc.θ)
#     Plots.plot(θ, z_centers)
# end
# Plots.mp4(anim, joinpath(path, "KM_θ.mp4"), fps = 10)

# anim = Plots.@animate for u in solver.u
#     qv = parent(u.Yc.qv)
#     Plots.plot(qv, z_centers)
# end
# Plots.mp4(anim, joinpath(path, "KM_qv.mp4"), fps = 10)

θ_end = parent(solver.u[end].Yc.θ)
qv_end = parent(solver.u[end].Yc.qv)
Plots.png(
    Plots.plot(θ_end, z_centers),
    joinpath(path, "KM_θ_end.png"),
)
Plots.png(
    Plots.plot(qv_end, z_centers),
    joinpath(path, "KM_qv_end.png"),
)
