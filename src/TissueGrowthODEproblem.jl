# Mechanical Relaxation ODE problem
function _fnc(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt = p
    for i = 1:N
        if i == 1
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,N],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,N],kₛ,l₀), τ(u[:,i+1],u[:,N]))*τ(u[:,i+1],u[:,N]) +
            Vₙ(u[:,N],u[:,i],u[:,i+1],kf,δt)
        elseif i == N
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,1],u[:,i-1],kₛ,l₀), τ(u[:,1],u[:,i-1]))*τ(u[:,1],u[:,i-1]) +
            Vₙ(u[:,i-1],u[:,i],u[:,1],kf,δt)
        else
            @views du[:,i] = (1/η) * dot(Fₛ⁺(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀) + Fₛ⁻(u[:,i],u[:,i+1],u[:,i-1],kₛ,l₀), τ(u[:,i+1],u[:,i-1]))*τ(u[:,i+1],u[:,i-1]) +
            Vₙ(u[:,i-1],u[:,i],u[:,i+1],kf,δt)
        end 
    end
end

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

# Normal Velocity ODE problem
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
function SetupODEproblem(btype,N,R₀,kₛ,η,kf,l₀,δt,Tmax)
    kₛ = kₛ*N
    η = η/N
    kf = kf/N
    # setting up initial conditions
    θ = collect(LinRange(0.0, 2*π, N+1));    # just use collect(θ) to convert into a vector
    pop!(θ);
    u0 = ElasticArray{Float64}(undef,2,N)
    for i in 1:N
        if btype == "circle"
            R = R₀ # to produce identical areas
            @views u0[:,i] .= [X(R,θ[i]), Y(R,θ[i])];
        elseif btype == "triangle"
            R = √((2*π*R₀^2)/sin(π/3))
            @views u0[:,i] .= [Xₜ(R,θ[i]*3/(2*π)), Yₜ(R,θ[i]*3/(2*π))];
        elseif btype == "square"
            R = √(π*(R₀^2)) # to produce identical areas
            @views u0[:,i] .= [Xₛ(R,θ[i]*2/pi), Yₛ(R,θ[i]*2/pi)];
        elseif btype == "hex"
            R = √((2/3*√3)*π*(R₀^2)) # to produce identical areas
            @views u0[:,i] .= [Xₕ(R,θ[i]*3/pi), Yₕ(R,θ[i]*3/pi)];
        end
    end
    #plotInitialCondition(u0)
    # solving ODE problem
    p = (N,kₛ,η,kf,l₀,δt)
    tspan = (0.0,Tmax)
    #prob = ODEProblem(_fnc2,u0,tspan,p)
    #prob = SplitODEProblem(_fnc1,_fnc2,u0,tspan,p)
    #return SplitODEProblem(_fnc1,_fnc2,u0,tspan,p)
    return ODEProblem(_fnc,u0,tspan,p)
end