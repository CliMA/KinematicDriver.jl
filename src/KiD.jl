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

const CP = CLIMAParameters
const APS = CP.AbstractEarthParameterSet

# Instantiate CliMA Parameters
struct AEPS <: CLIMAParameters.AbstractEarthParameterSet end
params = AEPS()

include("initial_condition.jl")
include("tendency.jl")
