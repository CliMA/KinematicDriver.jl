module Kinematic1D

import OrdinaryDiffEq
import NCDatasets
import UnPack
import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

import ClimaCore
import CLIMAParameters
import Thermodynamics
import CloudMicrophysics

const CC = ClimaCore
const CP = CLIMAParameters
const TD = Thermodynamics
const CM = CloudMicrophysics
const NC = NCDatasets
const ODE = OrdinaryDiffEq

include("EquationTypes.jl")
include("Kid_model.jl")
include("TimeStepping.jl")
include("NetCDFIO.jl")
include("callbacks.jl")
include("initial_condition.jl")
include("tendency.jl")
include("InternalParameters.jl")

end
