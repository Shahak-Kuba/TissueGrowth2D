function sim1D()
    # setting up simulation parameters
    N = 80 # number of cells
    m = 2 # number of springs per cell
    M = m*N # total number of springs along the interface
    R₀ = 1  # shape radius
    #kₛ = 0.7
    kₛ_Array = m*[0.1, 0.5, 1]
    l₀ = 1e-3/m
    kf = 0.7/m
    η = 1/m
    Tmax = 30 # days
    δt = 0.001
    #btypes = ["circle", "triangle", "square", "hex"]
    btype = "SineWave"
    savetimes = LinRange(0, Tmax, 30)

    #sol_array = Array{ODESolution}(undef,length(btypes));
    results = Vector{SimResults_t}(undef, 0)
    # creating 

    for ii in eachindex(kₛ_Array)
        @views kₛ = kₛ_Array[ii]
        prob, p = SetupODEproblem1D(btype, M, R₀, kₛ, η, kf, l₀, δt, Tmax)
        @time sol = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
        #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
        push!(results, postSimulation1D(btype, sol, p))
        printInfo(ii,length(kₛ_Array),btype,M,kₛ,η,kf)
    end

    return results

end