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

const CMNe = CM.MicrophysicsNonEq
const CM0 = CM.Microphysics0M
const CM1 = CM.Microphysics1M

include("EquationTypes.jl")
include("KiD_model.jl")
include("TimeStepping.jl")
include("NetCDFIO.jl")
include("callbacks.jl")
include("pysdm_functions.jl")
include("initial_condition.jl")
include("tendency.jl")
include("InternalParameters.jl")

end
