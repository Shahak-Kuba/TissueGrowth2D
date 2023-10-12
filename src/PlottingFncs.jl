
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

function plotResults(u, var)
    #f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    #    resolution=(500, 500))
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(1000, 400))
    ga = f[1, 1] = GridLayout()
    #gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), xlabel="x", ylabel="y")
    gaxmain = Axis(ga[1, 1], limits=(0, 10*pi, -1.5, 8), aspect=DataAspect(), xlabel="x", ylabel="y")
    CRange = findMinMax(var)
    for i in eachindex(u)
        lines!(gaxmain, u[i][1, :], u[i][2, :], color=var[i], colorrange=CRange,
            colormap=:jet, linewidth=3)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=:jet,
        flipaxis=false, label="ρ [μm²]")
    return f
end


function plotAreaVStime(sols)
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(700, 500))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], title="Void area over time", xlabel="Time [Days]", ylabel="Ω [μm²]")
    lin1 = lines!(gaxmain, sols[1].t, sols[1].Ω, linewidth=4, linestyle=:solid)
    lin2 = lines!(gaxmain, sols[2].t, sols[2].Ω, linewidth=4, linestyle=:dash)
    lin3 = lines!(gaxmain, sols[3].t, sols[3].Ω, linewidth=4, linestyle=:dot)
    lin4 = lines!(gaxmain, sols[4].t, sols[4].Ω, linewidth=4, linestyle=:dashdot)
    Legend(f[1, 2], [lin1, lin2, lin3, lin4], ["circle", "triangle", "square", "hex"])
    return f
end

function plotKapVsVel(sol)
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(700, 500))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], title="Curvature vs velocity @ different days", xlabel="κ", ylabel="vₙ [μms⁻¹]")
    for ii in axes(sol.Κ,1)
        if sol.btype == "triangle"
            mid = floor(Int64, length(sol.Κ[ii])/3)
        elseif sol.btype == "square"
            mid = floor(Int64, length(sol.Κ[ii])/4)
        elseif sol.btype == "hex"
            mid = floor(Int64, length(sol.Κ[ii])/6)
        elseif sol.btype == "circle"
            mid = floor(Int64, length(sol.Κ[ii])/2)
        end

        x = sort(sol.Κ[ii][2:mid])
        y = sort(sol.Vₙ[ii][2:mid])

        scatter!(gaxmain, x, y, markersize = 25, marker = '*')
        

        #fit = curve_fit(LogFit, x, y)
        #yfit = fit.(x)
        #lines!(gaxmain, x, yfit, linewidth=3, linestyle=:dash)
    end

    return f
end

