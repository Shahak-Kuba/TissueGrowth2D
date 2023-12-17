"""
SimResults_t(btype,t,u,∑F, ρ, V, Ω)

"""

struct SimResults_t
    btype::String
    t::Vector{Float64}
    u::Vector{Matrix{Float64}}
    ∑F::Vector{Vector{Float64}}
    Density::Vector{Matrix{Float64}}
    Vₙ::Vector{Vector{Float64}}
    Ω::Vector{Float64}
    ψ::Vector{Matrix{Float64}}
    Κ::Vector{Matrix{Float64}}
end