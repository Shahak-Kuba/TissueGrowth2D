# Mechanical Relaxation ODE problem
# Open boundary ODE Function (PERIODIC Boundary)
function ODE_fnc_1D_init!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt = p
    dom = 2*pi;
    @views for i in 1:N
        if i == 1
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N]-[dom, 0] ,kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N]-[dom, 0],kₛ,l₀), τ(u[:,i+1],u[:,N]-[dom, 0]))*τ(u[:,i+1],u[:,N]-[dom, 0])
        elseif i == N
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,1]+[dom, 0],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1]+[dom, 0],u[:,i-1],kₛ,l₀), τ(u[:,1]+[dom, 0],u[:,i-1]))*τ(u[:,1]+[dom, 0],u[:,i-1])
        else
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) 
        end 
    end
    nothing
end

"""

function ODE_fnc_1D!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt = p
    dom = 2*pi;
    @views for i in 1:N
        if i == 1
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N]-[dom, 0] ,kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N]-[dom, 0],kₛ,l₀), τ(u[:,i+1],u[:,N]-[dom, 0]))*τ(u[:,i+1],u[:,N]-[dom, 0]) +
            Vₙ(u[:,N]-[dom, 0],u[:,i],u[:,i+1],kf,δt)
        elseif i == N
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,1]+[dom, 0],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1]+[dom, 0],u[:,i-1],kₛ,l₀), τ(u[:,1]+[dom, 0],u[:,i-1]))*τ(u[:,1]+[dom, 0],u[:,i-1]) +
            Vₙ(u[:,i-1],u[:,i],u[:,1]+[dom, 0],kf,δt)
        else
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) +
            Vₙ(u[:,i-1],u[:,i],u[:,i+1],kf,δt)
        end 
    end
    nothing
end
"""

function ODE_fnc_1D!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt = p
    dom = 2*pi;
    uᵢ₊₁ = circshift(u,1)
    uᵢ₋₁ = circshift(u,-1)
    du[1,:] .= ((1/η) .* dot(Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁.-[dom 0],kₛ,l₀) + Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁.-[dom 0],kₛ,l₀), τ( uᵢ₊₁,uᵢ₋₁.-[dom 0]))*τ( uᵢ₊₁,uᵢ₋₁.-[dom 0]) +
                Vₙ(uᵢ₋₁.-[dom 0],u, uᵢ₊₁,kf,δt))[1,:]
    du[2:end-1,:] .= ((1/η) .* dot(Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀), τ(uᵢ₊₁,uᵢ₋₁))*τ(uᵢ₊₁,uᵢ₋₁) +
                Vₙ(uᵢ₋₁,u,uᵢ₊₁,kf,δt))[2:end-1,:]
    du[end,:] .= ((1/η) .* dot(Fₛ⁺(u,uᵢ₊₁.+[dom 0],uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u,uᵢ₊₁.+[dom 0],uᵢ₋₁,kₛ,l₀), τ(uᵢ₊₁.+[dom 0],uᵢ₋₁))*τ(uᵢ₊₁.+[dom 0],uᵢ₋₁) +
                Vₙ(uᵢ₋₁,u,uᵢ₊₁.+[dom 0],kf,δt))[end,:]
    nothing
end


function ODE_fnc_2D_init!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt = p
    @views for i in 1:N
        if i == 1
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N],kₛ,l₀), τ(u[:,i+1],u[:,N]))*τ(u[:,i+1],u[:,N])
        elseif i == N
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1],u[:,i-1],kₛ,l₀), τ(u[:,1],u[:,i-1]))*τ(u[:,1],u[:,i-1])
        else
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1])
        end 
    end
    nothing
end

"""
function ODE_fnc_2D!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt = p
    @views for i in 1:N
        if i == 1
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N],kₛ,l₀), τ(u[:,i+1],u[:,N]))*τ(u[:,i+1],u[:,N]) +
            Vₙ(u[:,N],u[:,i],u[:,i+1],kf,δt)
        elseif i == N
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1],u[:,i-1],kₛ,l₀), τ(u[:,1],u[:,i-1]))*τ(u[:,1],u[:,i-1]) +
            Vₙ(u[:,i-1],u[:,i],u[:,1],kf,δt)
        else
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) +
            Vₙ(u[:,i-1],u[:,i],u[:,i+1],kf,δt)
        end 
    end
    nothing
end
"""

function ODE_fnc_2D!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt = p
    uᵢ₊₁ = circshift(u,1)
    uᵢ₋₁ = circshift(u,-1)
    du .= (1/η) .* row_dot(Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀), τ(uᵢ₊₁,uᵢ₋₁)).*τ(uᵢ₊₁,uᵢ₋₁) +
                        Vₙ(uᵢ₋₁,u,uᵢ₊₁,kf,δt)
    """
    @views for i in 1:N
        if i == 1
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N],kₛ,l₀), τ(u[:,i+1],u[:,N]))*τ(u[:,i+1],u[:,N]) +
            Vₙ(u[:,N],u[:,i],u[:,i+1],kf,δt)
        elseif i == N
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1],u[:,i-1],kₛ,l₀), τ(u[:,1],u[:,i-1]))*τ(u[:,1],u[:,i-1]) +
            Vₙ(u[:,i-1],u[:,i],u[:,1],kf,δt)
        else
            @inbounds du[:,i] .= (1/η) .* dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) +
            Vₙ(u[:,i-1],u[:,i],u[:,i+1],kf,δt)
        end 
    end
    """
    nothing
end


"""
SetupODEproblem(btype,N,R₀,kₛ,η,kf,l₀,δt,Tmax) is a function to create the split ODE problem, for DifferentalEquations.jl. 

The inputs to this equations include:

btype: boundary type, 

N: number of springs in the system, 

R₀: radius for a circle based problem (this will get adjusted for different shapes), 

kₛ,η: the prescaled mechanical relaxation coefficients, 

kf: the amount of tissue produced per cell per unit time, 

l₀: resting spring length, 

δt: Euler timestep size, 

Tmax: end of simulation time
"""
function SetupODEproblem1D(btype,N,R₀,kₛ,η,kf,l₀,δt,Tmax)
    kₛ = kₛ*N
    η = η/N
    kf = kf/N
    # setting up initial conditions
    u0 = u0SetUp(btype,R₀,N)
    # solving ODE problem
    p = (N,kₛ,η,kf,l₀,δt)
    tspan = (0.0,Tmax)
    return ODEProblem(ODE_fnc_1D!,u0,tspan,p), p
end

function SetupODEproblem2D(btype,N,R₀,kₛ,η,kf,l₀,δt,Tmax)
    kₛ = kₛ*N
    η = η/N
    kf = kf/N
    u0 = u0SetUp(btype,R₀,N)
    #plotInitialCondition(u0)
    # solving ODE problem
    p = (N,kₛ,η,kf,l₀,δt)
    tspan = (0.0,Tmax)
    return ODEProblem(ODE_fnc_2D!,u0,tspan,p), p
end