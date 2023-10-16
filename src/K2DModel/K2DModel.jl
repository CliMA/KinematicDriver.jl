module K2DModel

import OrdinaryDiffEq
import NCDatasets

import ClimaCore
import Thermodynamics

const CC = ClimaCore
const TD = Thermodynamics
const NC = NCDatasets
const ODE = OrdinaryDiffEq

import ..K1DModel as K1D

include("./ode_utils.jl")
include("./netcdf_io.jl")
include("./initial_condition.jl")
include("./tendency.jl")

end
