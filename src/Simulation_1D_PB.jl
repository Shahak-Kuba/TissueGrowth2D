function sim1D_PB()
    # setting up simulation parameters
    #N = 384 # number of cells
    N = 250
    R₀ = 1  # shape radius
    #kₛ = 0.7
    kₛ_Array = [0.1, 0.7, 1]
    l₀ = 1e-3
    kf = 0.7
    η = 1
    Tmax = 50 # days
    δt = 0.001
    #btypes = ["circle", "triangle", "square", "hex"]
    btype = "SineWave"
    savetimes = LinRange(0, Tmax, 30)

    #sol_array = Array{ODESolution}(undef,length(btypes));
    results = Vector{SimResults_t}(undef, 0)
    # creating 

    for ii in eachindex(kₛ_Array)
        @views kₛ = kₛ_Array[ii]
        prob, p = SetupODEproblem1D_PB(btype, N, R₀, kₛ, η, kf, l₀, δt, Tmax)
        @time sol = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
        #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
        push!(results, postSimulation(btype, sol, p))
        printInfo(ii,length(kₛ_Array),btype,N,kₛ,η,kf)
    end

    return results

end