using Test

using ClimaCorePlots, Plots

include("../src/KiD.jl")

const FT = Float64

# Instantiate CliMA Parameters
struct AEPS <: CP.AbstractEarthParameterSet end
params = AEPS()

# Set up the computational domain and time step
z_min = FT(0)
z_max = FT(2e3)
n_elem = 128
Δt = 10.0
t_ini = 0.0
t_end = 10.0

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
w = Geometry.WVector.(ones(FT, face_space))

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

z_centers = parent(Fields.coordinate_field(space))
θ_liq_ice_end = parent(aux.θ_liq_ice)
q_tot_end = parent(solver.u[end].q_tot)

# Placeholder for testing
@testset "Test NaNs" begin
    @test !any(isnan.(θ_liq_ice_end))
    @test !any(isnan.(q_tot_end))
end
@testset "Test positive definite" begin
    @test minimum(θ_liq_ice_end) > FT(0)
    @test minimum(q_tot_end) > FT(0)
end

# Save some plots
ENV["GKSwstype"] = "nul"
Plots.GRBackend()

dir = "advect"
path = joinpath(@__DIR__, "output", dir)
mkpath(path)

anim = Plots.@animate for u in solver.u
    q_tot = parent(u.q_tot)
    Plots.plot(q_tot, z_centers)
end
Plots.mp4(anim, joinpath(path, "KM_qt.mp4"), fps = 10)

Plots.png(
    Plots.plot(θ_liq_ice_end, z_centers),
    joinpath(path, "KM_θ_end.png"),
)
Plots.png(
    Plots.plot(q_tot_end, z_centers),
    joinpath(path, "KM_qt_end.png"),
)
