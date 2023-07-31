function plotInitialCondition(u0)
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
        resolution = (1000, 700))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], xlabel = "x", ylabel = "y")
    append!(u0, u0[:,1])
    lines!(gaxmain,u0[1,:], u0[2,:], markersize = 10)
    return f
end

function plotResults(sol)
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
        resolution = (1000, 700))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-2,2,-2,2), aspect = DataAspect(), xlabel = "x", ylabel = "y")
    for i = 1:size(sol.u,1)
        append!(sol.u[i], sol.u[i][:,1])
        lines!(gaxmain,sol.u[i][1,:], sol.u[i][2,:], color = :blue, linewidth = 4)
    end
    return f
end

