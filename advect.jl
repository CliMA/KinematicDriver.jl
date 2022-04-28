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
z_max = FT(3e3)
n_elem = 128

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

# Thermodynamic states (TODO)
# - solve for density and pressure as a function of z 
function rhod(T,p,qv)
    return Thermodynamics.air_density(params, T, p, qv)
end

function rhod(θ,qv,rhod0)
    R_d::FT = CLIMAParameters.Planet.R_d(params)
    R_v::FT = CLIMAParameters.Planet.R_v(params)
    c_pv::FT = CLIMAParameters.Planet.cp_v(params)
    c_pd::FT = CLIMAParameters.Planet.cp_d(params)
    g::FT = CLIMAParameters.Planet.grav(params)
    
    function drho_dz(rhod, _, z)
        dql_dz=0
        θ_dry = θ * (1 + qv)^(R_d / c_pd) # TODO: fix this with thermodynamics.jl
        T = θ_dry * ((rhod * θ_dry) / 1000 * R_d)^(R_d/c_pd / (1 - R_d/c_pd))
        p = Thermodynamics.air_pressure_given_θ(params, T, rhod, qv) # TODO: should be moist rhod
        lv = Thermodynamics.latent_heat_vapor(params, T)
        Rq = R_v / (1/qv + 1) + R_d / (1+qv)
        cp = c_pv / (1/qv + 1) + c_pd / (1 + qv)
        rho = Thermodynamics.air_density(params, T, p)

        drhodz = (g / T * rho * (Rq / cp - 1) - p * lv / cp / T^2 * dql_dz) / Rq
        return drhodz
    end

    zspan = (z_min, z_max)
    prob = ODEProblem(drhodz, rhod0, (z_min, z_max))
    sol = solve(problem);

    return z -> sol.u(z)
end

# Define initial condition
function init_1d_column(::Type{FT}, params, z) where {FT}
    # physics parameters


    z_0::FT = 0.0
    z_1::FT = 740.0
    z_2::FT = 3260.0
    qv_0::FT = 0.016
    qv_1::FT = 0.0138
    qv_2::FT = 0.0024
    θ_0::FT = 279.9
    θ_1::FT = 279.9
    θ_2::FT = 312.66

    p0::FT = 1002 # hPa

    # potential temperature
    θ = z < z_1 ? θ_0 : θ_1 + (θ_2 - θ_1)/(z_2 - z_1) * z

    # water vapour specific humidity (TODO - or is it mixing ratio?)
    qv = z < z_1 ? qv_0 + (qv_1 - qv_0)/(z_1 - z_0) * z : qv_1 + (qv_2 - qv_1)/(z_2 - z_1) * z
    
    # density, pressure, etc...
    rhod0 = rhod(θ_0, p0, qv_0)
    rho = rhod(θ, p0, qv)

    return(θ = θ, qv = qv, rho=rho)
end

Yc = map(coord -> init_1d_column(FT, params, coord.z), coord)
w = Geometry.WVector.(zeros(FT, face_space))
Y = Fields.FieldVector(Yc = Yc, w = w)



# Advection Equation: ∂ϕ/dt = -∂(vΦ)
function advection_tendency!(dY, Y, _, t, w_params)
    Yc = Y.Yc

    w = Y.w
    w1 = w_params.w1
    t1 = w_params.t1
    w = t < t1 ? w1 * sin.(π * t / t1) : 0

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
          bottom = Operators.Extrapolate(),   #TODO - boundary conditions!
          top = Operators.Extrapolate(),
    )

    A_qv = Operators.AdvectionC2C(
           bottom = Operators.Extrapolate(),   #TODO - boundary conditions!
           top = Operators.Extrapolate(),
    )

    #@. dθ = -A_θ(w, θ) + fcc(w, θ)
    @. dθ = 0.0
    @. dqv = -A_qv(w, qv) + fcc(w, qv)
    return dY
end

# TODO - add more tendencies
# function microphysics_tendency!(dθ, θ, _, t)

# Solve the ODE operator
Δt = 10
problem = ODEProblem(advection_tendency!, Y, (0.0, 3600.0))
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
