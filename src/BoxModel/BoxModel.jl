module BoxModel

import CloudMicrophysics as CM

const CM1 = CM.Microphysics1M
const CM2 = CM.Microphysics2M
const CMP = CM.Parameters

include("parameters.jl")
import .Parameters as BP

include("equation_types.jl")
include("ode_utils.jl")
include("helper_functions.jl")
include("tendency.jl")

end
