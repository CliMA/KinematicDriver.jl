import ClimaCore:
    Fields,
    Domains,
    Topologies,
    Meshes,
    DataLayouts,
    Operators,
    Geometry,
    Spaces

import SurfaceFluxes
import Thermodynamics
import CloudMicrophysics
import CLIMAParameters

using OrdinaryDiffEq: ODEProblem, solve, SSPRK33

import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())
const FT = Float64

a = FT(0.0)
b = FT(4pi)
n = 128
α = FT(0.1)

w1 = 2.0
#t1 = 600.0 for future time varying w

domain = Domains.IntervalDomain(
    Geometry.ZPoint{FT}(a),
    Geometry.ZPoint{FT}(b),
    boundary_names = (:left, :right),
)
mesh = Meshes.IntervalMesh(domain, nelems = n)

cs = Spaces.CenterFiniteDifferenceSpace(mesh)
fs = Spaces.FaceFiniteDifferenceSpace(cs)

# VELOCITY: constant right now
V = w1 .* Geometry.WVector.(ones(FT, fs))

# Interpolate the water vapor and potTemp profiles
z_init = [0.0, 740.0, 3260.0]
qv_init = [0.016, 0.0138, 0.0024]
theta_init = [279.9, 279.9, 312.66]

# interpolate between the I.C. points
function interpolate_z(z_init, q_init, z_interp)
    # find the two spanning z points
    for i=1:length(z_init)-1
        if z_interp > z_init[i] && z_interp < z_init[i+1]
            z1 = z_init[i]
            i1 = i
            z2 = z_init[i+1]
            i2 = i+1
            break
        end
        z1, z2 .= z_init[end]
        i1, i2 .= length(z_init)
    end

    # do the linear interpolation
    q1, q2 = q_init[i1], q_init[i2]
    if z1 != z2
        q_interp = q1 + (q2 - q1)/(z2 - z1)*(z_interp - z1)
    else
        q_interp = q1
    end

    return q_interp
end

print(Fields.coordinate_field(cs).z)
qv_initz = z -> interpolate_z(z_init, qv_init, z) where {z <: FT}
qv = qv_initz.(Fields.coordinate_field(cs).z)

θ_initz = z -> interpolate_z(z_init, theta_init, z) where {z <: FT}
θ = θ_initz.(Fields.coordinate_field(cs).z)


# TRACER(S) that gets advected: theta and qt (up for negotiation)
θ = sin.(Fields.coordinate_field(cs).z)
qt = sin.(Fields.coordinate_field(cs).z)


# Solve advection Equation: ∂θ/dt = -∂(vθ)
# use the advection operator
function tendency3!(dθ, θ, _, t)

    fcc = Operators.FluxCorrectionC2C(
        left = Operators.Extrapolate(),
        right = Operators.Extrapolate(),
    )
    fcf = Operators.FluxCorrectionF2F(
        left = Operators.Extrapolate(),
        right = Operators.Extrapolate(),
    )
    A = Operators.AdvectionC2C(
        left = Operators.SetValue(sin(-t)),
        right = Operators.Extrapolate(),
    )
    return @. dθ = -A(V, θ)
end



# # Solve the ODE operator
# Δt = 0.001
# prob3 = ODEProblem(tendency3!, θ, (0.0, 10.0))
# sol3 = solve(
#     prob3,
#     SSPRK33(),
#     dt = Δt,
#     saveat = 10 * Δt,
#     progress = true,
#     progress_message = (dt, u, p, t) -> t,
# );

# ENV["GKSwstype"] = "nul"
# using ClimaCorePlots, Plots
# Plots.GRBackend()

# dir = "advect"
# path = joinpath(@__DIR__, "output", dir)
# mkpath(path)


# anim = Plots.@animate for u in sol3.u
#     Plots.plot(u, xlim = (-1, 1))
# end
# Plots.mp4(anim, joinpath(path, "C2C_advect.mp4"), fps = 10)
# Plots.png(
#     Plots.plot(sol3.u[end], xlim = (-1, 1)),
#     joinpath(path, "sol3_advect_end.png"),
# )

# function linkfig(figpath, alt = "")
#     # buildkite-agent upload figpath
#     # link figure in logs if we are running on CI
#     if get(ENV, "BUILDKITE", "") == "true"
#         artifact_url = "artifact://$figpath"
#         print("\033]1338;url='$(artifact_url)';alt='$(alt)'\a\n")
#     end
# end

# linkfig(
#     relpath(joinpath(path, "advect_end.png"), joinpath(@__DIR__, "../..")),
#     "Advect End Simulation",
# )
