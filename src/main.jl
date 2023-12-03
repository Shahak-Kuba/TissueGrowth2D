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
sols1D = sim1D();
f = plotResults1D(sols1D[1].u, sols1D[1].Density)
f = plotResults1D(sols1D[2].u, sols1D[2].Density)
f = plotResults1D(sols1D[3].u, sols1D[3].Density)
f = plotResults1D(sols1D[4].u, sols1D[4].Density)

f = plotResults1D_Velocity(sols1D[1].u, sols1D[1].Vₙ)
f = plotResults1D_Velocity(sols1D[2].u, sols1D[2].Vₙ)
f = plotResults1D_Velocity(sols1D[3].u, sols1D[3].Vₙ)
f = plotResults1D_Velocity(sols1D[4].u, sols1D[4].Vₙ)
#f = plotKapVsVel(sols[1])


# 2D simulation 
sols2D = sim2D();
f = plotAreaVStime(sols2D)

cmap = :viridis
f = plotResults2D(sols2D[1][1].u, sols2D[1][1].Density, cmap)
f = plotResults2D(sols2D[1][2].u, sols2D[1][2].Density, cmap)
f = plotResults2D(sols2D[1][3].u, sols2D[1][3].Density, cmap)
f = plotResults2D(sols2D[1][4].u, sols2D[1][4].Density, cmap)
f = plotResults2D(sols2D[1][5].u, sols2D[1][5].Density, cmap)
f = plotResults2D(sols2D[1][6].u, sols2D[1][6].Density, cmap)

f = plotResults2D(sols2D[2][1].u, sols2D[2][1].Density, cmap)
f = plotResults2D(sols2D[2][2].u, sols2D[2][2].Density, cmap)
f = plotResults2D(sols2D[2][3].u, sols2D[2][3].Density, cmap)
f = plotResults2D(sols2D[2][4].u, sols2D[2][4].Density, cmap)
f = plotResults2D(sols2D[2][5].u, sols2D[2][5].Density, cmap)
f = plotResults2D(sols2D[2][6].u, sols2D[2][6].Density, cmap)

f = plotResults2D(sols2D[3][1].u, sols2D[3][1].Density, cmap)
f = plotResults2D(sols2D[3][2].u, sols2D[3][2].Density, cmap)
f = plotResults2D(sols2D[3][3].u, sols2D[3][3].Density, cmap)
f = plotResults2D(sols2D[3][4].u, sols2D[3][4].Density, cmap)
f = plotResults2D(sols2D[3][5].u, sols2D[3][5].Density, cmap)
f = plotResults2D(sols2D[3][6].u, sols2D[3][6].Density, cmap)

f = plotKapVsVel(sols2D[1])
f = plotKapVsVel(sols2D[2])
f = plotKapVsVel(sols2D[3])
f = plotKapVsVel(sols2D[4])
f = plotKapVsVel(sols2D[5])
