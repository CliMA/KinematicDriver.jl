module Kinematic1D

import ClimaCore
import Thermodynamics
import CloudMicrophysics
import CLIMAParameters

import NCDatasets
import OrdinaryDiffEq
import UnPack

import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

const CC = ClimaCore
const TD = Thermodynamics
const CP = CLIMAParameters
const NC = NCDatasets
const ODE = OrdinaryDiffEq

include("EquationTypes.jl")
include("Kid_model.jl")
include("TimeStepping.jl")
include("NetCDFIO.jl")
include("callbacks.jl")
include("initial_condition.jl")
include("tendency.jl")

end
