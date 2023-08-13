using OrdinaryDiffEq
using CairoMakie
using ColorSchemes
using Colors
using ElasticArrays
using LinearAlgebra
using BenchmarkTools
using Printf

include("MechanicalEqns.jl")
include("PoreBoundaries.jl")
include("PlottingFncs.jl")
include("GeometrySolvers.jl")
include("TissueGrowthODEproblem.jl")
include("Misc.jl")
include("DataStructs.jl")

function main2()
    # setting up simulation parameters
    N = 384 # number of cells
    R₀ = 1  # shape radius
    kₛ = 0.025   # high Fₛ: 2.5, mid Fₛ: 0.5, low Fₛ: 0.01 
    l₀ = 1e-3
    kf = 5e-2
    η = 1
    Tmax = 60 # days
    δt = 0.001
    btypes = ["circle", "triangle", "square", "hex"]
    savetimes = LinRange(0, Tmax, 8)

    #sol_array = Array{ODESolution}(undef,length(btypes));
    results = Vector{SimResults_t}(undef, 0)
    # creating 

    for ii in eachindex(btypes)
        @views btype = btypes[ii]
        prob, p = SetupODEproblem(btype, N, R₀, kₛ, η, kf, l₀, δt, Tmax)
        #@time sol = solve(prob,SplitEuler(),saveat=savetimes,dt=δt)
        @time sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
        push!(results, postSimulation(btype, sol, p))
        #printInfo(ii,length(btypes),btype,N,kₛ,η,kf)
    end

    return results

end

sols = main2();
f = plotAreaVStime(sols)

f = plotResults(sols[1].u, sols[1].ψ)
f = plotResults(sols[2].u, sols[2].ψ)
f = plotResults(sols[3].u, sols[3].ψ)
f = plotResults(sols[4].u, sols[4].ψ)

f = plotKapVsVel(sols[1])
f = plotKapVsVel(sols[2])
f = plotKapVsVel(sols[3])
f = plotKapVsVel(sols[4])
