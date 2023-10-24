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
include("PlottingFncs.jl")
include("GeometrySolvers.jl")
include("TissueGrowthODEproblem.jl")
include("Misc.jl")
include("DataStructs.jl")
include("Simulation_1D.jl")
include("Simulation_1D_PB.jl")
include("Simulation_2D.jl")


"""
# 1D simulation (non periodic boundary)
sols = sim1D();
#f = plotAreaVStime(sols)
f = plotResults(sols[1].u, sols[1].Density)
f = plotKapVsVel(sols[1])



# 1D Simulation (free right boundary)
sols = sim1D_FB();
#f = plotAreaVStime(sols)
f = plotResults(sols[1].u, sols[1].Density)
f = plotKapVsVel(sols[1])
"""


# 1D simulation (periodic boundary)
sols = sim1D_PB();
f = plotResults(sols[1].u, sols[1].Density)
f = plotResults(sols[2].u, sols[2].Density)
f = plotResults(sols[3].u, sols[3].Density)
f = plotResults(sols[4].u, sols[4].Density)
#f = plotKapVsVel(sols[1])


# 2D simulation 
sols = sim2D();
f = plotAreaVStime(sols)

f = plotResults(sols[1].u, sols[1].Density)
f = plotResults(sols[2].u, sols[2].Density)
f = plotResults(sols[3].u, sols[3].Density)
f = plotResults(sols[4].u, sols[4].Density)

f = plotKapVsVel(sols[1])
f = plotKapVsVel(sols[2])
f = plotKapVsVel(sols[3])
f = plotKapVsVel(sols[4])
