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

Vₙ(ρ⁺,ρ⁻,kf) = kf*(ρ⁺+ρ⁻)/2  

"""
ϕ is the vector projection scalar
"""

ϕ(u,v) = (dot(u,v)/(norm(v)^2))

"""
"""
function ∠uv(rₗ, rₘ, rᵣ) 
   u = rᵣ - rₘ
   v = rₗ - rₘ
   ψ = round(dot(u,v)/(norm(u)*norm(v)), digits=2)
   return acos(ψ)
end

"""
ξ is the scaling function to fix the danger zone problem cause by ∠ < π/2 
function takes inputs left node: rₗ (i-1), central node: rₘ (i), and right node: rᵣ (i+1)
"""
function ξ(rₗ, rₘ, rᵣ)
    @views ∠ = ∠uv(rₗ, rₘ, rᵣ);
    # linear mapping of scaling values ξ = m∠ + c
    #m = (1-√2)/(π/2)
    #c = 1 - m*π
    return (1-√2)/(π/2)*∠ + 1 - ((1-√2)/(π/2))*π
end