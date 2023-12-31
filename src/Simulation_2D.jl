function sim2D()
    # setting up simulation parameters
    N = 384/4 # number of cells
    m = 2 # number of springs per cell
    M = Int(m*N) # total number of springs along the interface
    #N = 500
    R₀ = 1  # shape radius
    #ks_Array = [0.01, 0.25, 1] # 0.25 is nice smoothing
    D = [0.001]#, 0.075, 0.15, 1]
    l₀ = 1
    kf = 0.001
    η = 1
    growth_dir = "inward" 
    Tmax = 28# days
    δt = 0.0005
    btypes = ["star", "cross"]#["circle", "triangle", "square", "hex", "star","cross"]
    savetimes = LinRange(0, Tmax, 8)

    all_results = Vector{Vector{SimResults_t}}(undef, 0)

    for jj in eachindex(D)
        @views kₛ = D[jj]*(η)/((l₀)^2)
        #sol_array = Array{ODESolution}(undef,length(btypes));
        results = Vector{SimResults_t}(undef, 0)
        # creating 

        for ii in eachindex(btypes)
            @views btype = btypes[ii]
            prob, p = SetupODEproblem2D(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir)
            @time sol = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
            #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
            push!(results, postSimulation2D(btype, sol, p))
            printInfo(ii,length(btypes),btype,N,kₛ*m,η/m,kf/m,M,D[jj])
        end
        push!(all_results,results)
    end

    return all_results

end


function sim2D_δt()
    # setting up simulation parameters
    N = 384/4 # number of cells
    m = 2 # number of springs per cell
    M = Int(m*N) # total number of springs along the interface
    #N = 500
    R₀ = 1  # shape radius
    kₛ  = 0.1 # 0.25 is nice smoothing
    l₀ = 1e-3
    kf = 0.1
    η = 1
    Tmax = 28# days
    δt_Array = [0.00001,0.0001,0.001,0.005]
     #["circle", "triangle", "square", "hex", "star","cross"]
    btypes = ["triangle", "square", "hex", "star","cross"]
    savetimes = LinRange(0, Tmax, 8)

    all_results = Vector{Vector{SimResults_t}}(undef, 0)

    for jj in eachindex(btypes)
        results = Vector{SimResults_t}(undef, 0)
        @views btype = btypes[jj]

        for ii in eachindex(δt_Array)
            @views δt = δt_Array[ii]
            # creating 
            prob, p = SetupODEproblem2D(btype, M, R₀, kₛ, η, kf, l₀, δt, Tmax)
            @time sol = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
            #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
            push!(results, postSimulation2D(btype, sol, p))
            printInfo(ii,length(btypes),btype,M,kₛ*M,η/M,kf/M)
        end

        push!(all_results,results)
    end

    return all_results

end