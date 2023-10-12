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
function main()
    # setting up simulation parameters
    #N = 384 # number of cells
    N = 500
    R₀ = 1  # shape radius
    kₛ = 0.01
    l₀ = 1e-3
    kf = 1
    η = 1
    Tmax = 250 # days
    δt = 0.001
    #btypes = ["circle", "triangle", "square", "hex"]
    btypes = ["SineWave"]
    savetimes = LinRange(0, Tmax, 30)

    #sol_array = Array{ODESolution}(undef,length(btypes));
    results = Vector{SimResults_t}(undef, 0)
    # creating 

    for ii in eachindex(btypes)
        @views btype = btypes[ii]
        prob, p = SetupODEproblem(btype, N, R₀, kₛ, η, kf, l₀, δt, Tmax)
        @time sol = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
        #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
        push!(results, postSimulation(btype, sol, p))
        printInfo(ii,length(btypes),btype,N,kₛ,η,kf)
    end

    return results

end
"""

# 1D simulation (non periodic boundary)
sols = sim1D();
#f = plotAreaVStime(sols)
f = plotResults(sols[1].u, sols[1].Density)
f = plotKapVsVel(sols[1])



# 1D Simulation (free right boundary)
sols = sim1D_FB();
f = plotAreaVStime(sols)
f = plotResults(sols[1].u, sols[1].Density)
f = plotKapVsVel(sols[1])


# 1D simulation (periodic boundary)
sols = sim1D_PB();
f = plotAreaVStime(sols)

f = plotResults(sols[1].u, sols[1].Density)
f = plotResults(sols[2].u, sols[2].Density)
f = plotResults(sols[3].u, sols[3].Density)
f = plotResults(sols[4].u, sols[4].Density)

f = plotKapVsVel(sols[1])
f = plotKapVsVel(sols[2])
f = plotKapVsVel(sols[3])
f = plotKapVsVel(sols[4])


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
