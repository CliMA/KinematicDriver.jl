module K2DModel

import NCDatasets as NC

import ClimaCore as CC
import Thermodynamics as TD

import ..Common as CO
import ..K1DModel as K1D

include("./ode_utils.jl")
include("./initial_condition.jl")
include("./tendency.jl")

end
