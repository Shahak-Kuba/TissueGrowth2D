using OrdinaryDiffEq
using CairoMakie
using ElasticArrays
using LinearAlgebra
using BenchmarkTools

include("DataStructs.jl")
include("GeometrySolvers.jl")
include("MechanicalEqns.jl")
include("Statistics.jl")
include("Callbacks.jl")
include("PoreBoundaries.jl")
include("PlottingFncs.jl")
include("Misc.jl")
include("TissueGrowthODEproblem.jl")

function main()
    # setting up simulation parameters
    N = 120 # number of cells
    R = 10  # shape radius
    kₛ = 0.1   # high Fₛ: 2.5, mid Fₛ: 0.5, low Fₛ: 0.01 
    l₀ = 1e-3
    kf = 1e-4
    η = 1


    # rescaling 
    kₛ = kₛ*N
    η = η/N
    kf = kf/N

    Tmax = 1 # days
    btype = "hex"

    # setting up initial conditions
    θ = collect(LinRange(0.0, 2*π, N+1));    # just use collect(θ) to convert into a vector
    pop!(θ);
    u0 = ElasticArray{Float64}(undef,2,N)
    for i in 1:N
        if btype == "circle"
            @views u0[:,i] .= [X(R,θ[i]), Y(R,θ[i])];
        elseif btype == "square"
            @views u0[:,i] .= [Xₛ(R,θ[i]*2/pi), Yₛ(R,θ[i]*2/pi)];
        elseif btype == "hex"
            @views u0[:,i] .= [Xₕ(R,θ[i]*3/pi), Yₕ(R,θ[i]*3/pi)];
        end
    end
    #plotInitialCondition(u0)
    # solving ODE problem
    p = (N,kₛ,η,kf,l₀)
    tspan = (0.0,Tmax)
    prob = ODEProblem(_fnc2,u0,tspan,p)
    savetimes = LinRange(0,Tmax,10)

    u0, sol = solve(prob,BS5(),saveat=savetimes)
    # plotting

    return sol
end

u0, sol = @time main();

f = plotResults(sol)