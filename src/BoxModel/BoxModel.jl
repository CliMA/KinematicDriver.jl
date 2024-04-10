module BoxModel

import CloudMicrophysics as CM

const CM1 = CM.Microphysics1M
const CM2 = CM.Microphysics2M
const CMP = CM.Parameters

import ..Common as CO

include("parameters.jl")
import .Parameters as BP

include("ode_utils.jl")
include("tendency.jl")

end
