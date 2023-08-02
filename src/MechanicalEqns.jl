"""
δ is the length between two nodes in the positive direction
"""
δ(rᵢ₊₁, rᵢ) = √((rᵢ₊₁[1]-rᵢ[1])^2 + (rᵢ₊₁[2]-rᵢ[2])^2)

"""
τ calculates the unit tangent vector between two neighbouring points rᵢ₊₁ and rᵢ₋₁ of a central point rᵢ
"""
τ(rᵢ₊₁, rᵢ₋₁) = (rᵢ₊₁ - rᵢ₋₁)/δ(rᵢ₊₁, rᵢ₋₁)

"""
n calculates the unit normal vector between two neighbouring points rᵢ₊₁ and rᵢ₋₁ of a central point rᵢ
"""
n(rᵢ₊₁, rᵢ₋₁) = [-τ(rᵢ₊₁, rᵢ₋₁)[2];τ(rᵢ₊₁, rᵢ₋₁)[1]]


"""
Fₛ⁺ is the spring force (Nonlinear) used for mechanical relaxation in the positive direction.
l = length of the spring
l₀ = resting length of the spring
kₛ = spring coefficient
"""

Fₛ⁺(rᵢ,rᵢ₊₁, rᵢ₋₁,kₛ,l₀) = kₛ*l₀^2 * (1/l₀ - 1/δ(rᵢ₊₁, rᵢ)) * τ(rᵢ₊₁, rᵢ)

"""
Fₛ⁻ is the spring force (Nonlinear) used for mechanical relaxation in the negative direction.
l = length of the spring
l₀ = resting length of the spring
kₛ = spring coefficient
"""

Fₛ⁻(rᵢ,rᵢ₊₁, rᵢ₋₁,kₛ,l₀) = -kₛ*l₀^2 * (1/l₀ - 1/δ(rᵢ, rᵢ₋₁)) * τ(rᵢ, rᵢ₋₁)

"""
"""
ρ(rᵢ₊₁, rᵢ) = 1/δ(rᵢ₊₁, rᵢ);

"""
Vₙ is the normal velocity of the interface such that Vₙ ∝ ρ
ρ = the density of the cell
kf = the amount of tissue produced per unit area per unit time
"""

#Vₙ(ρ⁺::Float64,ρ⁻::Float64,kf::Float64) = kf*(ρ⁺+ρ⁻)/2  

function Vₙ(rᵢ₋₁,rᵢ,rᵢ₊₁,kf,δt)
    ρₗ = ρ(rᵢ,rᵢ₋₁)
    ρᵣ = ρ(rᵢ₊₁,rᵢ)
    Vₗ = kf*ρₗ
    Vᵣ = kf*ρᵣ

    nₗ = n(rᵢ₋₁,rᵢ)
    nᵣ = n(rᵢ,rᵢ₊₁)

    rₘ₁ = rᵢ - Vₗ*nₗ*δt
    rₗ = rᵢ₋₁ - Vₗ*nₗ*δt
    rₘ₂ = rᵢ - Vᵣ*nᵣ*δt
    rᵣ = rᵢ₊₁ - Vᵣ*nᵣ*δt

    return (lineIntersection(rₘ₁,rₗ,rₘ₂,rᵣ) - rᵢ)/δt
end


"""
κ(rᵢ₋₁,rᵢ,rᵢ₊₁) approximated the curvature of the shape for κ vs V plots
"""

function κ(rᵢ₋₁,rᵢ,rᵢ₊₁)
    
    @views Dt2X = 0.5*(rᵢ₋₁[1] -2*rᵢ[1] + rᵢ₊₁[1]);
    @views DtX = 0.5*(rᵢ₊₁[1] - rᵢ₋₁[1]);
    @views Dt2Y = 0.5*(rᵢ₋₁[2] -2*rᵢ[2] + rᵢ₊₁[2]);
    @views DtY = 0.5*(rᵢ₊₁[2] - rᵢ₋₁[2]);

    return (DtX*Dt2Y - DtY*Dt2X) / (DtX^2 + DtY^2)^(3/2); 
end