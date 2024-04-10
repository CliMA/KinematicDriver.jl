module K1DModel

import OrdinaryDiffEq as ODE
import NCDatasets as NC
import UnPack
import Logging
import TerminalLoggers

import ClimaCore as CC
import Thermodynamics as TD
import CloudMicrophysics as CM

const CMP = CM.Parameters
const CMNe = CM.MicrophysicsNonEq
const CM0 = CM.Microphysics0M
const CM1 = CM.Microphysics1M
const CM2 = CM.Microphysics2M
const CMAM = CM.AerosolModel
const CMAA = CM.AerosolActivation

import ..Common as CO

include("parameters.jl")
import .Parameters as KP

include("ode_utils.jl")
include("netcdf_io.jl")
include("callbacks.jl")
include("helper_functions.jl")
include("initial_condition.jl")
include("tendency.jl")

Logging.global_logger(TerminalLoggers.TerminalLogger())

end
