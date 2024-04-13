module Common

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

include("parameters.jl")
include("equation_types.jl")
include("helper_functions.jl")
include("ode_utils.jl")
include("pysdm_functions.jl")
include("tendency.jl")

end