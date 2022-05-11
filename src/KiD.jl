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

import NCDatasets
import OrdinaryDiffEq
import UnPack

import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

const TD = Thermodynamics
const CP = CLIMAParameters
const NC = NCDatasets
const ODE = OrdinaryDiffEq

# Instantiate CliMA Parameters
struct AEPS <: CP.AbstractEarthParameterSet end
params = AEPS()

include("TimeStepping.jl")
include("NetCDFIO.jl")
include("callbacks.jl")
include("initial_condition.jl")
include("tendency.jl")
