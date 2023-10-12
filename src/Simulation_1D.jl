function sim1D()
    # setting up simulation parameters
    #N = 384 # number of cells
    N = 300
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

    #sol_array = Array{ODESolution}(undef,length(bty1pes));
    results = Vector{SimResults_t}(undef, 0)
    # creating 

    for ii in eachindex(btypes)
        @views btype = btypes[ii]
        prob, p = SetupODEproblem1D(btype, N, R₀, kₛ, η, kf, l₀, δt, Tmax)
        @time sol = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
        #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
        push!(results, postSimulation(btype, sol, p))
        printInfo(ii,length(btypes),btype,N,kₛ,η,kf)
    end

    return results

end