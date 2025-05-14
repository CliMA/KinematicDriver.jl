module Common

import OrdinaryDiffEq as ODE
import SpecialFunctions as SF
import NCDatasets as NC
using Statistics

import CSV as CSV
import Interpolations as IP

import ClimaCore as CC
import Thermodynamics as TD
import CloudMicrophysics as CM
import Cloudy as CL

const CMP = CM.Parameters
const CMNe = CM.MicrophysicsNonEq
const CM0 = CM.Microphysics0M
const CM1 = CM.Microphysics1M
const CM2 = CM.Microphysics2M
const CMD = CM.CloudDiagnostics
const CMAM = CM.AerosolModel
const CMAA = CM.AerosolActivation
const CMP3 = CM.P3Scheme

include("parameters.jl")
include("equation_types.jl")
include("helper_functions.jl")
include("ode_utils.jl")
include("netcdf_io.jl")
include("initial_condition.jl")
include("pysdm_functions.jl")
include("tendency.jl")

end
