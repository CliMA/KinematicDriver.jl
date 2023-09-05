module Kinematic1D

import OrdinaryDiffEq as ODE
import NCDatasets as NC
import UnPack
import Logging
import TerminalLoggers

import ClimaCore as CC
import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics:
    CommonTypes as CMT,
    MicrophysicsNonEq as CMNe,
    Microphysics0M as CM0,
    Microphysics1M as CM1,
    Microphysics2M as CM2,
    AerosolModel as CMAM,
    AerosolActivation as CMAA

# modules used only for calibrateKiD
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

# include driveKiD and calibrateKiD
include("driveKiD/driveKiD.jl")
include("calibrateKiD/calibrateKiD.jl")
Logging.global_logger(TerminalLoggers.TerminalLogger())

end
