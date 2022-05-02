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

using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, Tsit5

import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

const TD = Thermodynamics
const CP = CLIMAParameters

# Instantiate CliMA Parameters
struct AEPS <: CP.AbstractEarthParameterSet end
params = AEPS()

include("initial_condition.jl")
include("tendency.jl")
