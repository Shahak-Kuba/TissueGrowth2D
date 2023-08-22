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

function main()
    # setting up simulation parameters
    N = 96 # number of cells
    R₀ = 1  # shape radius
    kₛ = 1  
    l₀ = 1e-3
    kf = 1.4e-1
    η = 1
    Tmax = 22 # days
    δt = 0.001
    btypes = ["circle", "triangle", "square", "hex"]
    savetimes = LinRange(0, Tmax, 7)

    #sol_array = Array{ODESolution}(undef,length(btypes));
    results = Vector{SimResults_t}(undef, 0)
    # creating 

    for ii in eachindex(btypes)
        @views btype = btypes[ii]
        prob, p = SetupODEproblem(btype, N, R₀, kₛ, η, kf, l₀, δt, Tmax)
        #@time sol = solve(prob,SplitEuler(),saveat=savetimes,dt=δt)
        @time sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
        push!(results, postSimulation(btype, sol, p))
        printInfo(ii,length(btypes),btype,N,kₛ,η,kf)
    end

    return results

end

sols = main();
f = plotAreaVStime(sols)

f = plotResults(sols[1].u, sols[1].Density)
f = plotResults(sols[2].u, sols[2].Density)
f = plotResults(sols[3].u, sols[3].Density)
f = plotResults(sols[4].u, sols[4].Density)

f = plotKapVsVel(sols[1])
f = plotKapVsVel(sols[2])
f = plotKapVsVel(sols[3])
f = plotKapVsVel(sols[4])
