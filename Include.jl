# setup paths -
const _ROOT = pwd();
const _PATH_TO_SRC = joinpath(_ROOT, "src");
const _PATH_TO_DATA = joinpath(_ROOT, "data");
const _PATH_TO_CONFIG = joinpath(_ROOT, "config");
const _PATH_TO_NETWORK = joinpath(_ROOT, "network");
const _PATH_TO_SIMULATIONS = joinpath(_ROOT, "simulations");
const _PATH_TO_FIGS = joinpath(_ROOT, "figs");

# load external packages -
using DataFrames
using CSV
using Statistics
using LinearAlgebra
using NumericalIntegration
using Plots
using Colors
using DataInterpolations
using Interpolations
using JSON
using DifferentialEquations
using GlobalSensitivity
using ProgressMeter
using Symbolics
using DelimitedFiles

# load my codes -
include(joinpath(_PATH_TO_SRC, "Types.jl"))
include(joinpath(_PATH_TO_SRC, "Data.jl"))
include(joinpath(_PATH_TO_SRC, "Balances.jl"))
include(joinpath(_PATH_TO_SRC, "Kinetics.jl"))
include(joinpath(_PATH_TO_SRC, "Control.jl"))
include(joinpath(_PATH_TO_SRC, "Error.jl"))
include(joinpath(_PATH_TO_SRC, "Utility.jl"))
include(joinpath(_PATH_TO_SRC, "SolveBalances.jl"))
