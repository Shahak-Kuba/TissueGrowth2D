"""
function _fnc2(du,u,p,t) 
    N,kₛ,η,kf,l₀ = p
    for i = 1:N
        if i == 1
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N],kₛ,l₀), τ(u[:,i+1],u[:,N]))*τ(u[:,i+1],u[:,N]) + 
            Vₙ(u[:,N-1],u[:,N],u[:,i],u[:,i+1],u[:,i+2],kf)
        elseif i == 2
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) + 
            Vₙ(u[:,N],u[:,i-1],u[:,i],u[:,i+1],u[:,i+2],kf)
        elseif i == N-1
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) + 
            Vₙ(u[:,i-2],u[:,i-1],u[:,i],u[:,i+1],u[:,1],kf)
        elseif i == N
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1],u[:,i-1],kₛ,l₀), τ(u[:,1],u[:,i-1]))*τ(u[:,1],u[:,i-1]) + 
            Vₙ(u[:,i-2],u[:,i-1],u[:,i],u[:,1],u[:,2],kf)
        else
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) + 
            Vₙ(u[:,i-2],u[:,i-1],u[:,i],u[:,i+1],u[:,i+2],kf)
        end 
    end
end"""

function _fnc1(du,u,p,t) 
    N,kₛ,η,kf,l₀ = p
    for i = 1:N
        if i == 1
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N],kₛ,l₀), τ(u[:,i+1],u[:,N]))*τ(u[:,i+1],u[:,N])
        elseif i == N
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1],u[:,i-1],kₛ,l₀), τ(u[:,1],u[:,i-1]))*τ(u[:,1],u[:,i-1]) 
        else
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) 
        end 
    end
end

function _fnc2(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt = p
    for i = 1:N
        if i == 1
            @views du[:,i] = Vₙ(u[:,N],u[:,i],u[:,i+1],kf,δt)
        elseif i == N
            @views du[:,i] = Vₙ(u[:,i-1],u[:,i],u[:,1],kf,δt)
        else
            @views du[:,i] = Vₙ(u[:,i-1],u[:,i],u[:,i+1],kf,δt)
        end 
    end
end