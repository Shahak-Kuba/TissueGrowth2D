using OrdinaryDiffEq
using CairoMakie
using ColorSchemes
using Colors
using ElasticArrays
using LinearAlgebra
using BenchmarkTools
using Printf
using CurveFit

include("MechanicalEqns.jl")
include("PoreBoundaries.jl")
include("PlottingFncs1D.jl")
include("PlottingFncs2D.jl")
include("Simulation_1D.jl")
include("Simulation_2D.jl")
include("TissueGrowthODEproblem.jl")
include("GeometrySolvers.jl")
include("Misc.jl")
include("DataStructs.jl")
include("PostSimulation.jl")



# 1D simulation (periodic boundary)
sols = sim1D();
f = plotResults1D(sols[1].u, sols[1].Density)
f = plotResults1D(sols[2].u, sols[2].Density)
f = plotResults1D(sols[3].u, sols[3].Density)
f = plotResults1D(sols[4].u, sols[4].Density)

f = plotResults1D_Velocity(sols[1].u, sols[1].Vₙ)
f = plotResults1D_Velocity(sols[2].u, sols[2].Vₙ)
f = plotResults1D_Velocity(sols[3].u, sols[3].Vₙ)
f = plotResults1D_Velocity(sols[4].u, sols[4].Vₙ)
#f = plotKapVsVel(sols[1])


# 2D simulation 
sols = sim2D();
f = plotAreaVStime(sols)

cmap = :viridis
f = plotResults2D(sols[1].u, sols[1].Density, cmap)
f = plotResults2D(sols[2].u, sols[2].Density, cmap)
f = plotResults2D(sols[3].u, sols[3].Density, cmap)
f = plotResults2D(sols[4].u, sols[4].Density, cmap)

f = plotKapVsVel(sols[1])
f = plotKapVsVel(sols[2])
f = plotKapVsVel(sols[3])
f = plotKapVsVel(sols[4])
