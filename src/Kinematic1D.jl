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
const CMAM = CM.AerosolModel
const CMAA = CM.AerosolActivation

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
