
function printInfo(simNum,simTotal,btype,N,kₛ,η,kf)
    println(@sprintf "----------------------------- Simulation %d/%d Complete -----------------------------" simNum simTotal)
    println(@sprintf "Boundary Type: %s, Cell count: %d, kₛ: %.5f, η: %.5f, kf: %.5f" btype N kₛ η kf)
    println(@sprintf "-----------------------------------------------------------------------------------")
end