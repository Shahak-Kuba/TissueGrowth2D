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

# setting up ODE problem
function _fnc(du,u,p,t) 
    N,kₛ,η,kf,l₀ = p
    for i = 1:N
        if i == 1
            du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N],kₛ,l₀), τ(u[:,i+1],u[:,N]))*τ(u[:,i+1],u[:,N]) + Vₙ(ρ(u[:,i+1],u[:,i]), ρ(u[:,i], u[:,N]), kf)*n(u[:,i+1],u[:,N])
        elseif i == N
            du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1],u[:,i-1],kₛ,l₀), τ(u[:,1],u[:,i-1]))*τ(u[:,1],u[:,i-1]) + Vₙ(ρ(u[:,1],u[:,i]), ρ(u[:,i], u[:,i-1]), kf)*n(u[:,1],u[:,i-1])
        else
            du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) + Vₙ(ρ(u[:,i+1],u[:,i]), ρ(u[:,i], u[:,i-1]), kf)*n(u[:,i+1],u[:,i-1])
        end 
    end
end

function main()

    # setting up simulation parameters
    N = 100 # number of cells
    R = 0.1  # shape radius
    kₛ = 0.05
    l₀ = 1e-3
    kf = 1.9e-3;
    η = 1;
    # setting up initial conditions
    θ = collect(LinRange(0.0, 2*π, N+1));    # just use collect(θ) to convert into a vector
    pop!(θ);
    u0 = ElasticArray{Float64}(undef,2,N)
    for i in 1:N
        @views u0[:,i] .= [X(R,θ[i]), Y(R,θ[i])];
    end
    # solving ODE problem
    p = (N,kₛ,η,kf,l₀)
    tspan = (0.0,40.0)
    prob = ODEProblem(_fnc,u0,tspan,p)
    savetimes = LinRange(1,40,8)

    sol = solve(prob,Tsit5(),saveat=savetimes)
    # plotting

    return sol
end

sol = @time main();

f = plotResults(sol)
