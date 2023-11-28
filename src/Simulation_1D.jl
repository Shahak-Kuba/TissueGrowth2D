function sim1D()
    # setting up simulation parameters
    N = 163 # number of cells
    m = 1 # number of springs per cell
    M = m*N # total number of springs along the interface
    R₀ = 1  # shape radius
    kₛ_Array = [13]
    #kₛ_Array = [0.1, 1, 2]
    l₀ = 1e-3
    kf = 0.7
    η = 1
    Tmax = 50 # days
    δt = 0.0001
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