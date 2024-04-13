module K1DModel

import OrdinaryDiffEq as ODE
import SpecialFunctions as SF
import Logging
import TerminalLoggers

import ClimaCore as CC
import Thermodynamics as TD
import CloudMicrophysics as CM

const CMAM = CM.AerosolModel
const CMAA = CM.AerosolActivation

import ..Common as CO

include("parameters.jl")
include("ode_utils.jl")
include("initial_condition.jl")
include("tendency.jl")

Logging.global_logger(TerminalLoggers.TerminalLogger())

end
