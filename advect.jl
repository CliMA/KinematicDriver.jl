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

# Instantiate CliMA Parameters
struct AEPS <: CLIMAParameters.AbstractEarthParameterSet end
params = AEPS()

# Set up the computational domain
z_min = FT(0)
z_max = FT(2e3)
n_elem = 128

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

# Define initial condition
function init_1d_column(::Type{FT}, params, z) where {FT}
    # physics parameters
    # R_d::FT = CLIMAParameters.Planet.R_d(params)

    z_0::FT = 0.0
    z_1::FT = 740.0
    z_2::FT = 3260.0
    qv_0::FT = 0.016
    qv_1::FT = 0.0138
    qv_2::FT = 0.0024
    θ_0::FT = 279.9
    θ_1::FT = 279.9
    θ_2::FT = 312.66

    # density, pressure, etc...

    # potential temperature
    θ = z < z_1 ? θ_0 : θ_1 + (θ_2 - θ_1)/(z_2 - z_1) * z

    # water vapour specific humidity (TODO - or is it mixing ratio?)
    qv = z < z_1 ? qv_0 + (qv_1 - qv_0)/(z_1 - z_0) * z : qv_1 + (qv_2 - qv_1)/(z_2 - z_1) * z

    return(θ = θ, qv = qv)
end

Yc = map(coord -> init_1d_column(FT, params, coord.z), coord)
w = Geometry.WVector.(ones(FT, face_space)) #TODO - should be changing in time
Y = Fields.FieldVector(Yc = Yc, w = w)

# Advection Equation: ∂ϕ/dt = -∂(vΦ)
function advection_tendency!(dY, Y, _, t)

    Yc = Y.Yc
    w = Y.w
    # TODO @. w = Y.w * sin(t) ?

    dYc = dY.Yc

    θ = Yc.θ
    qv = Yc.qv

    dθ = dYc.θ
    dqv = dYc.qv

    fcc = Operators.FluxCorrectionC2C(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    fcf = Operators.FluxCorrectionF2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    A_θ = Operators.AdvectionC2C(
          bottom = Operators.SetValue(279.9),   #TODO - boundary conditions!
          top = Operators.Extrapolate(),
    )

    A_qv = Operators.AdvectionC2C(
           bottom = Operators.SetValue(0.016),   #TODO - boundary conditions!
           top = Operators.Extrapolate(),
    )

    @. dθ = -A_θ(w, θ)
    @. dqv = -A_qv(w, qv)
    return dY
end

# TODO - add more tendencies
# function microphysics_tendency!(dθ, θ, _, t)

# Solve the ODE operator
Δt = 0.001
problem = ODEProblem(advection_tendency!, Y, (0.0, 10.0))
solver = solve(
    problem,
    SSPRK33(),
    dt = Δt,
    saveat = 10 * Δt,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);

ENV["GKSwstype"] = "nul"
using ClimaCorePlots, Plots
Plots.GRBackend()

dir = "advect"
path = joinpath(@__DIR__, "output", dir)
mkpath(path)

z_centers = parent(Fields.coordinate_field(space))

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
