module CalibrateCMP

import OrdinaryDiffEq as ODE
import NCDatasets as NC
import ClimaCore as CC

using Distributions
using LinearAlgebra
using Random
using Plots
using Interpolations
using Optim
using EnsembleKalmanProcesses
using EnsembleKalmanProcesses.ParameterDistributions
using EnsembleKalmanProcesses.Observations
using JLD2

const EKP = EnsembleKalmanProcesses
const DS = Distributions
const FT = Float64

import ..Common as CO
import ..BoxModel as BX
import ..K1DModel as KD

include("ReferenceStats.jl")
include("ReferenceModels.jl")
include("DistributionUtils.jl")
include("KiDUtils.jl")
include("OptimizationUtils.jl")
include("IOUtils.jl")
include("HelperFuncs.jl")
include("Diagnostics.jl")

end
