"""
δ is the length between two nodes in the positive direction
"""
δ(rᵢ₊₁, rᵢ) = √((rᵢ₊₁[1] - rᵢ[1])^2 + (rᵢ₊₁[2] - rᵢ[2])^2)

"""
τ calculates the unit tangent vector between two neighbouring points rᵢ₊₁ and rᵢ₋₁ of a central point rᵢ
"""
τ(rᵢ₊₁, rᵢ₋₁) = (rᵢ₊₁ .- rᵢ₋₁) / δ(rᵢ₊₁, rᵢ₋₁)

"""
n calculates the unit normal vector between two neighbouring points rᵢ₊₁ and rᵢ₋₁ of a central point rᵢ
"""
n(rᵢ₊₁, rᵢ₋₁) = [-τ(rᵢ₊₁, rᵢ₋₁)[2]; τ(rᵢ₊₁, rᵢ₋₁)[1]]
#n(rᵢ₊₁, rᵢ₋₁) = (τv = τ(rᵢ₊₁, rᵢ₋₁);return (-τv[2], τv[1]))

"""
Fₛ⁺ is the spring force (Nonlinear) used for mechanical relaxation in the positive direction.
l = length of the spring
l₀ = resting length of the spring
kₛ = spring coefficient
"""

Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) = kₛ * l₀^2 * (1 / l₀ - 1 / δ(rᵢ₊₁, rᵢ)) * τ(rᵢ₊₁, rᵢ)

"""
Fₛ⁻ is the spring force (Nonlinear) used for mechanical relaxation in the negative direction.
l = length of the spring
l₀ = resting length of the spring
kₛ = spring coefficient
"""

Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) = -kₛ * l₀^2 * (1 / l₀ - 1 / δ(rᵢ, rᵢ₋₁)) * τ(rᵢ, rᵢ₋₁)

function PostCalcs(u, p)
    N, kₛ, η, kf, l₀, δt = p

    ∑F = zeros(size(u, 2))
    density = zeros(size(u, 2))
    vₙ = zeros(size(u, 2))
    ψ = zeros(size(u, 2))
    Κ = zeros(size(u, 2))

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

    return ∑F, vₙ, density, ψ, Κ
end

"""

"""
ρ(rᵢ₊₁, rᵢ) = 1 / δ(rᵢ₊₁, rᵢ);

"""
Vₙ is the normal velocity of the interface such that Vₙ ∝ ρ
ρ = the density of the cell
kf = the amount of tissue produced per unit area per unit time
"""

#Vₙ(ρ⁺::Float64,ρ⁻::Float64,kf::Float64) = kf*(ρ⁺+ρ⁻)/2  

function Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt)
    ρₗ = ρ(rᵢ, rᵢ₋₁)
    ρᵣ = ρ(rᵢ₊₁, rᵢ)
    Vₗ = kf * ρₗ
    Vᵣ = kf * ρᵣ

    nₗ = n(rᵢ₋₁, rᵢ)
    nᵣ = n(rᵢ, rᵢ₊₁)

    rₘ₁ = rᵢ - Vₗ * nₗ * δt
    rₗ = rᵢ₋₁ - Vₗ * nₗ * δt
    rₘ₂ = rᵢ - Vᵣ * nᵣ * δt
    rᵣ = rᵢ₊₁ - Vᵣ * nᵣ * δt

    return (lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ) - rᵢ) / δt
end


"""
κ(rᵢ₋₁,rᵢ,rᵢ₊₁) approximated the curvature of the shape for κ vs V plots
"""

function κ(rᵢ₋₁, rᵢ, rᵢ₊₁)

    @views Dt2X = 0.5 * (rᵢ₋₁[1] - 2 * rᵢ[1] + rᵢ₊₁[1])
    @views DtX = 0.5 * (rᵢ₊₁[1] - rᵢ₋₁[1])
    @views Dt2Y = 0.5 * (rᵢ₋₁[2] - 2 * rᵢ[2] + rᵢ₊₁[2])
    @views DtY = 0.5 * (rᵢ₊₁[2] - rᵢ₋₁[2])

    return abs((DtX * Dt2Y - DtY * Dt2X) / (DtX^2 + DtY^2)^(3 / 2))
end

""" 
postSimulation()

Function to perform post simulation calculations which returns a data structure which contains all data
"""

function postSimulation(btype, sol, p)

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
        Fnet, nV, den, stre, kap = PostCalcs(sol.u[ii], p)
        push!(∑F, Fnet)
        push!(vₙ, nV)
        push!(density, den)
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_t(btype, sol.t, sol.u, ∑F, density, vₙ, Area, ψ, Κ)
end