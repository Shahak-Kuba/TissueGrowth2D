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

function main()
    # setting up simulation parameters
    N = 384 # number of cells
    R₀ = 1  # shape radius
    kₛ = 0.025   # high Fₛ: 2.5, mid Fₛ: 0.5, low Fₛ: 0.01 
    l₀ = 1e-3
    kf = 5e-2
    η = 1


    # rescaling 
    kₛ = kₛ*N
    η = η/N
    kf = kf/N

    Tmax = 60 # days
    δt = 0.001;
    btype = "triangle"

    # setting up initial conditions
    θ = collect(LinRange(0.0, 2*π, N+1));    # just use collect(θ) to convert into a vector
    pop!(θ);
    u0 = ElasticArray{Float64}(undef,2,N)
    for i in 1:N
        if btype == "circle"
            R = R₀ # to produce identical areas
            @views u0[:,i] .= [X(R,θ[i]), Y(R,θ[i])];
        elseif btype == "triangle"
            R = √((2*π*R₀^2)/sin(π/3))  # side length to produce identical areas
            @views u0[:,i] .= [Xₜ(R,θ[i]*3/(2*π)), Yₜ(R,θ[i]*3/(2*π))];
        elseif btype == "square"
            R = √(π*(R₀^2))  # side length to produce identical areas
            @views u0[:,i] .= [Xₛ(R,θ[i]*2/pi), Yₛ(R,θ[i]*2/pi)];
        elseif btype == "hex"
            R = √((2/(3√3))*π*(R₀^2))  # side length to produce identical areas
            @views u0[:,i] .= [Xₕ(R,θ[i]*3/pi), Yₕ(R,θ[i]*3/pi)];
        end
    end
    #plotInitialCondition(u0)
    # solving ODE problem
    p = (N,kₛ,η,kf,l₀,δt)
    tspan = (0.0,Tmax)
    #prob = ODEProblem(_fnc2,u0,tspan,p)
    prob = SplitODEProblem(_fnc1,_fnc2,u0,tspan,p)
    savetimes = LinRange(0,Tmax,8)

    sol = solve(prob,SplitEuler(),saveat=savetimes,dt=δt)
    # plotting

    return sol
end

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
    btypes = ["circle","triangle","square","hex"]
    savetimes = LinRange(0,Tmax,8)

    #sol_array = Array{ODESolution}(undef,length(btypes));
    results = Vector{SimResults_t}(undef,0)
    # creating 

    for ii in eachindex(btypes)
        @views btype = btypes[ii]
        prob, p = SetupODEproblem(btype,N,R₀,kₛ,η,kf,l₀,δt,Tmax)
        #@time sol = solve(prob,SplitEuler(),saveat=savetimes,dt=δt)
        @time sol = solve(prob,Euler(),saveat=savetimes,dt=δt)
        push!(results, postSimulation(btype, sol, p))
        #printInfo(ii,length(btypes),btype,N,kₛ,η,kf)
    end

    return results

end

sols = main2();
f = plotAreaVStime(sols)

u0, sol = @time main();

results = Vector{SimResults_t}(undef,0);
push!(results, postSimulation(btype, sol, p));

f = plotResults(sols[4])