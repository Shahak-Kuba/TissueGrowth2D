function PostCalcs1D(u, p)
    N, kₛ, η, kf, l₀, δt = p

    ∑F = zeros(size(u, 1))
    density = zeros(size(u, 1))
    vₙ = zeros(size(u, 1))
    ψ = zeros(size(u, 1))
    Κ = zeros(size(u, 1))

    """

    for i in axes(u, 2)
        if i == 1
            @views ∑F[i] = abs((1 / η) * dot(Fₛ⁺(u[:, i], u[:, i+1], u[:, N]-[2*pi, 0], kₛ, l₀) + Fₛ⁻(u[:, i], u[:, i+1], u[:, N]-[2*pi, 0], kₛ, l₀), τ(u[:, i+1], u[:, N]-[2*pi, 0])))
            @views vₙ[i] = norm(Vₙ(u[:,N]-[2*pi, 0],u[:,i],u[:,i+1],kf,δt))
            @views density[i] = ρ(u[:, N]-[2*pi, 0], u[:, i])
            @views ψ[i] = ∑F[i] / δ(u[:, i+1], u[:, i])
            @views Κ[i] = κ(u[:, N], u[:, i], u[:, i+1]-[2*pi, 0])
        elseif i == N
            @views ∑F[i] = abs((1 / η) * dot(Fₛ⁺(u[:, i], u[:, 1]+[2*pi, 0], u[:, i-1], kₛ, l₀) + Fₛ⁻(u[:, i], u[:, 1]+[2*pi, 0], u[:, i-1], kₛ, l₀), τ(u[:, 1]+[2*pi, 0], u[:, i-1])))
            @views vₙ[i] = norm(Vₙ(u[:,i-1],u[:,i],u[:,1]+[2*pi, 0],kf,δt))
            @views density[i] = ρ(u[:, 1]+[2*pi, 0], u[:, i])
            @views ψ[i] = ∑F[i] / δ(u[:, 1], u[:, i])
            @views Κ[i] = κ(u[:, i-1], u[:, i], u[:, 1]+[2*pi, 0])
        elseif i == N + 1
            continue
        else
            @views ∑F[i] = abs((1 / η) * dot(Fₛ⁺(u[:, i], u[:, i+1], u[:, i-1], kₛ, l₀) + Fₛ⁻(u[:, i], u[:, i+1], u[:, i-1], kₛ, l₀), τ(u[:, i+1], u[:, i-1])))
            @views vₙ[i] = norm(Vₙ(u[:, i-1], u[:, i], u[:, i+1], kf, δt))
            @views density[i] = ρ(u[:, i+1], u[:, i])
            @views ψ[i] = ∑F[i] / δ(u[:, i+1], u[:, i])
            @views Κ[i] = κ(u[:, i-1], u[:, i], u[:, i+1])
        end
    end
    """

    uᵢ₊₁ = circshift(u,1)
    uᵢ₋₁ = circshift(u,-1)

    ∑F = dot(Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀), τ(uᵢ₊₁,uᵢ₋₁))


    return ∑F, vₙ, density, ψ, Κ
end


""" 
postSimulation1D()

Function to perform post simulation calculations which returns a data structure which contains all data
"""

function postSimulation1D(btype, sol, p)

    c = size(sol.t, 1)

    Area = Vector{Float64}(undef, c)
    ∑F = Vector{Vector{Float64}}(undef, 0)
    ψ = Vector{Vector{Float64}}(undef, 0)
    density = Vector{Vector{Float64}}(undef, 0)
    vₙ = Vector{Vector{Float64}}(undef, 0)
    Κ = Vector{Vector{Float64}}(undef, 0)

    for ii in axes(sol.u, 1)
        Area[ii] = Ω(sol.u[ii]) # area calculation
        #append!(sol.u[ii], sol.u[ii][:,1]) # closing the domain Ω
        Fnet, nV, den, stre, kap = PostCalcs1D(sol.u[ii], p)
        push!(∑F, Fnet)
        push!(vₙ, nV)
        push!(density, den)
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_t(btype, sol.t, sol.u, ∑F, density, vₙ, Area, ψ, Κ)
end


###################################################################################################


function PostCalcs2D(u, p)
    N, kₛ, η, kf, l₀, δt = p

    ∑F = zeros(size(u, 1))
    density = zeros(size(u, 1))
    vₙ = zeros(size(u, 1))
    ψ = zeros(size(u, 1))
    Κ = zeros(size(u, 1))

    """
    for i in axes(u, 2)
        if i == 1
            @views ∑F[i] = abs((1 / η) * dot(Fₛ⁺(u[:, i], u[:, i+1], u[:, N], kₛ, l₀) + Fₛ⁻(u[:, i], u[:, i+1], u[:, N], kₛ, l₀), τ(u[:, i+1], u[:, N])))
            @views vₙ[i] = norm(Vₙ(u[:, N], u[:, i], u[:, i+1], kf, δt))
            @views density[i] = ρ(u[:, i+1], u[:, i])
            @views ψ[i] = ∑F[i] / δ(u[:, i+1], u[:, i])
            @views Κ[i] = κ(u[:, N], u[:, i], u[:, i+1])
        elseif i == N
            @views ∑F[i] = abs((1 / η) * dot(Fₛ⁺(u[:, i], u[:, 1], u[:, i-1], kₛ, l₀) + Fₛ⁻(u[:, i], u[:, 1], u[:, i-1], kₛ, l₀), τ(u[:, 1], u[:, i-1])))
            @views vₙ[i] = norm(Vₙ(u[:, i-1], u[:, i], u[:, 1], kf, δt))
            @views density[i] = ρ(u[:, 1], u[:, i])
            @views ψ[i] = ∑F[i] / δ(u[:, 1], u[:, i])
            @views Κ[i] = κ(u[:, i-1], u[:, i], u[:, 1])
        elseif i == N + 1
            continue
        else
            @views ∑F[i] = abs((1 / η) * dot(Fₛ⁺(u[:, i], u[:, i+1], u[:, i-1], kₛ, l₀) + Fₛ⁻(u[:, i], u[:, i+1], u[:, i-1], kₛ, l₀), τ(u[:, i+1], u[:, i-1])))
            @views vₙ[i] = norm(Vₙ(u[:, i-1], u[:, i], u[:, i+1], kf, δt))
            @views density[i] = ρ(u[:, i+1], u[:, i])
            @views ψ[i] = ∑F[i] / δ(u[:, i+1], u[:, i])
            @views Κ[i] = κ(u[:, i-1], u[:, i], u[:, i+1])
        end
    end
    """

    uᵢ₊₁ = circshift(u,1)
    uᵢ₋₁ = circshift(u,-1)
    ∑F = row_dot(Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀), τ(uᵢ₊₁,uᵢ₋₁))
    density = ρ(uᵢ₊₁, u)
    ψ = ∑F ./ δ(uᵢ₊₁, u)
    #Κ[i] = κ(u[:, i-1], u[:, i], u[:, i+1])

    return ∑F, vₙ, density, ψ, Κ
end

""" 
postSimulation2D()

Function to perform post simulation calculations which returns a data structure which contains all data
"""

function postSimulation2D(btype, sol, p)

    c = size(sol.t, 1)

    Area = Vector{Float64}(undef, c)
    ∑F = Vector{Vector{Float64}}(undef, 0)
    ψ = Vector{Vector{Float64}}(undef, 0)
    density = Vector{Vector{Float64}}(undef, 0)
    vₙ = Vector{Vector{Float64}}(undef, 0)
    Κ = Vector{Vector{Float64}}(undef, 0)

    for ii in axes(sol.u, 1)
        Area[ii] = Ω(sol.u[ii]) # area calculation
        #append!(sol.u[ii], sol.u[ii][:,1]) # closing the domain Ω
        Fnet, nV, den, stre, kap = PostCalcs2D(sol.u[ii], p)
        push!(∑F, vcat(Fnet...))
        push!(vₙ, nV)
        push!(density, vcat(den...))
        push!(ψ, vcat(stre...))
        push!(Κ, kap)
    end

    return SimResults_t(btype, sol.t, sol.u, ∑F, density, vₙ, Area, ψ, Κ)
end