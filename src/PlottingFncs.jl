
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


# Colormaps available at: https://docs.juliahub.com/MakieGallery/Ql23q/0.2.17/generated/colors.html#Colormaps

function plotResults(u,var)
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
        resolution = (500, 500))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-1.5,1.5,-1.5,1.5), aspect = DataAspect(), xlabel = "x", ylabel = "y")
    CRange = findMinMax(var)
    for i in eachindex(u)
        lines!(gaxmain,u[i][1,:], u[i][2,:], color = var[i], colorrange = CRange,
        colormap = :jet, linewidth = 4)
    end
    Colorbar(f[1,2], limits = CRange, colormap = :jet,
    flipaxis = false, label = "ρ [μm²]")
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



