function plotResults(sol)
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
        resolution = (500, 500))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-2,2,-2,2), aspect = DataAspect(), xlabel = "x", ylabel = "y")
    for i = 1:size(sol.u,1)
        #append!(sol.u[i], sol.u[i][:,1])
        #append!(sol.ψ[i], sol.ψ[i][:,1])
        lines!(gaxmain,sol.u[i][1,:], sol.u[i][2,:], color = sol.Vₙ[i], colorrange = (0,0.01),
        colormap = :RdYlBu_7, linewidth = 4)
    end
    
    return f
end


function plotAreaVStime(sols)
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
        resolution = (500, 500))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], xlabel = "Time [Days]", ylabel = "Ω [μm²]")
    for ii in eachindex(sols)
        lines!(gaxmain, sols[ii].t, sols[ii].Ω, linewidth = 4)
    end
    return f
end



