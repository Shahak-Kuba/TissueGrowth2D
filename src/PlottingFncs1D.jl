# code to automatically save figures into a folder

function findMinMax(var)
    Min = 0
    Max = 0
    for ii in eachindex(var)
        if ii == 1
            Min = minimum(var[ii])
            Max = maximum(var[ii])
        else
            if minimum(var[ii]) < Min
                Min = minimum(var[ii])
            end
            if maximum(var[ii]) > Max
                Max = maximum(var[ii])
            end
        end
    end
    return (Min, Max)
end

function plotResults1D(u, var)
    #f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    #    resolution=(500, 500))
    f = Figure(fontsize = 32,backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    #gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), xlabel="x", ylabel="y")
    gaxmain = Axis(ga[1, 1], limits=(0, 2*pi, 1, 8), aspect=DataAspect(), xlabel="x", ylabel="y")
    CRange = findMinMax(var)
    #CRange = (0.05,0.2)
    for i in eachindex(u)
        lines!(gaxmain, u[i][:,1], u[i][:,2], color=var[i].data, colorrange=CRange,
            colormap=:jet, linewidth=5)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=:jet,
        flipaxis=false, label="v [mm/day]") 
    return f
end

function plotResults1D_Velocity(u, var)
    #f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    #    resolution=(500, 500))
    f = Figure(fontsize = 32,backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    #gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), xlabel="x", ylabel="y")
    gaxmain = Axis(ga[1, 1], limits=(0, 2*pi, 1, 8), aspect=DataAspect(), xlabel="x", ylabel="y")
    CRange = findMinMax(var)
    #CRange = (0.05,0.2)
    for i in eachindex(u)
        lines!(gaxmain, u[i][:,1], u[i][:,2], color=var[i], colorrange=CRange,
            colormap=:jet, linewidth=5)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=:jet,
        flipaxis=false, label="v [mm/day]") 
    return f
end