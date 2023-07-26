using OrdinaryDiffEq
using CairoMakie
using ElasticArrays
using LinearAlgebra
using BenchmarkTools

include("DataStructs.jl")
include("MechanicalEqns.jl")
include("Statistics.jl")
include("Callbacks.jl")
include("PoreBoundaries.jl")
include("PlottingFncs.jl")
include("Misc.jl")
include("TissueGrowthODEproblem.jl")

# setting up ODE problem
function _fnc(du,u,p,t) 
    N,kₛ,η,kf,l₀ = p
    for i = 1:N
        if i == 1
            du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N],kₛ,l₀), τ(u[:,i+1],u[:,N]))*τ(u[:,i+1],u[:,N]) + 
            ξ(u[:,i+1],u[:,i],u[:,N])*Vₙ(ρ(u[:,i+1],u[:,i]), ρ(u[:,i], u[:,N]), kf)*n(u[:,i+1],u[:,N])
        elseif i == N
            du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1],u[:,i-1],kₛ,l₀), τ(u[:,1],u[:,i-1]))*τ(u[:,1],u[:,i-1]) + 
            ξ(u[:,1],u[:,i],u[:,i-1])*Vₙ(ρ(u[:,1],u[:,i]), ρ(u[:,i], u[:,i-1]), kf)*n(u[:,1],u[:,i-1])
        else
            du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) + 
            ξ(u[:,i+1],u[:,i],u[:,i-1])*Vₙ(ρ(u[:,i+1],u[:,i]), ρ(u[:,i], u[:,i-1]), kf)*n(u[:,i+1],u[:,i-1])
        end 
    end
end

function main()
    # setting up simulation parameters
    N = 120 # number of cells
    R = 1  # shape radius
    kₛ = 1
    l₀ = 1e-3
    kf = 2.2e-4
    η = 1

    # rescaling 
    kₛ = kₛ*N
    η = η/N
    #kf = kf/N

    Tmax = 60 # days
    btype = "hex"

    # setting up initial conditions
    θ = collect(LinRange(0.0, 2*π, N+1));    # just use collect(θ) to convert into a vector
    pop!(θ);
    u0 = ElasticArray{Float64}(undef,2,N)
    for i in 1:N
        if btype == "circle"
            @views u0[:,i] .= [X(R,θ[i]), Y(R,θ[i])];
        elseif btype == "square"
            @views u0[:,i] .= [Xₛ(R,θ[i]*4/pi), Yₛ(R,θ[i]*4/pi)];
        elseif btype == "hex"
            @views u0[:,i] .= [Xₕ(R,θ[i]*3/pi), Yₕ(R,θ[i]*3/pi)];
        end
    end
    #plotInitialCondition(u0)
    # solving ODE problem
    p = (N,kₛ,η,kf,l₀)
    tspan = (0.0,Tmax)
    prob = ODEProblem(_fnc,u0,tspan,p)
    savetimes = LinRange(0,Tmax,7)

    u0, sol = solve(prob,BS5(),saveat=savetimes)
    # plotting

    return sol
end

u0, sol = @time main();

f = plotResults(sol)