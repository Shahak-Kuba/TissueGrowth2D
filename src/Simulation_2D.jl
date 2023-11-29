function sim2D()
    # setting up simulation parameters
    N = 384/4 # number of cells
    m = 2 # number of springs per cell
    M = Int(m*N) # total number of springs along the interface
    #N = 500
    R₀ = 1  # shape radius
    kₛ = 0.25
    l₀ = 1e-3
    kf = 0.1
    η = 1
    Tmax = 31# days
    δt = 0.0005
    btypes = ["circle", "triangle", "square", "hex", "star"]
    savetimes = LinRange(0, Tmax, 8)

    #sol_array = Array{ODESolution}(undef,length(btypes));
    results = Vector{SimResults_t}(undef, 0)
    # creating 

    for ii in eachindex(btypes)
        @views btype = btypes[ii]
        prob, p = SetupODEproblem2D(btype, M, R₀, kₛ, η, kf, l₀, δt, Tmax)
        @time sol = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
        #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
        push!(results, postSimulation2D(btype, sol, p))
        printInfo(ii,length(btypes),btype,M,kₛ*M,η/M,kf/M)
    end

    return results

end