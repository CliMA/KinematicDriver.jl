using Test

using ClimaCorePlots, Plots

include("../KiD.jl")

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
w = Geometry.WVector.(ones(FT, face_space)) #TODO - should be changing in time
Y = Fields.FieldVector(Yc = Yc, w = w)

# Solve the ODE operator
problem = ODEProblem(advection_tendency!, Y, (t_ini, t_end))
solver = solve(
    problem,
    SSPRK33(),
    dt = Δt,
    saveat = 10 * Δt,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);

z_centers = parent(Fields.coordinate_field(space))
θ_end = parent(solver.u[end].Yc.θ)
qv_end = parent(solver.u[end].Yc.qv)

# Placeholder for testing
@testset "Test NaNs" begin
    @test !any(isnan.(θ_end))
    @test !any(isnan.(qv_end))
end
@testset "Test positive definite" begin
    @test minimum(θ_end) > FT(0)
    @test minimum(qv_end) > FT(0)
end

# Save some plots
ENV["GKSwstype"] = "nul"
Plots.GRBackend()

dir = "advect"
path = joinpath(@__DIR__, "output", dir)
mkpath(path)

anim = Plots.@animate for u in solver.u
    θ = parent(u.Yc.θ)
    Plots.plot(θ, z_centers)
end
Plots.mp4(anim, joinpath(path, "KM_θ.mp4"), fps = 10)

anim = Plots.@animate for u in solver.u
    qv = parent(u.Yc.qv)
    Plots.plot(qv, z_centers)
end
Plots.mp4(anim, joinpath(path, "KM_qv.mp4"), fps = 10)

Plots.png(
    Plots.plot(θ_end, z_centers),
    joinpath(path, "KM_θ_end.png"),
)
Plots.png(
    Plots.plot(qv_end, z_centers),
    joinpath(path, "KM_qv_end.png"),
)
