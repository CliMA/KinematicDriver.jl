module Kinematic1D

import OrdinaryDiffEq
import NCDatasets
import UnPack
import Logging
import TerminalLoggers

import ClimaCore
import Thermodynamics
import CloudMicrophysics

const CC = ClimaCore
const TD = Thermodynamics
const CM = CloudMicrophysics
const NC = NCDatasets
const ODE = OrdinaryDiffEq

const CMT = CM.CommonTypes
const CMNe = CM.MicrophysicsNonEq
const CM0 = CM.Microphysics0M
const CM1 = CM.Microphysics1M
const CM2 = CM.Microphysics2M

include("Parameters.jl")
import .Parameters as KP

include("EquationTypes.jl")
include("KiD_model.jl")
include("TimeStepping.jl")
include("NetCDFIO.jl")
include("callbacks.jl")
include("helper_functions.jl")
include("dispatch_helper.jl")
include("pysdm_functions.jl")
include("initial_condition.jl")
include("tendency.jl")

include("calibration_pipeline/calibrateKiD.jl")
Logging.global_logger(TerminalLoggers.TerminalLogger())

end
