function plotResults(sol)
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
        resolution = (1000, 700))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], xlabel = "x", ylabel = "y")
    for i = 1:size(sol.u,1)
        append!(sol.u[i], sol.u[i][:,1])
        lines!(gaxmain,sol.u[i][1,:], sol.u[i][2,:])
    end
    return f
end