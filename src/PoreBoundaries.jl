
# Circular Boundary
X(R,θ) = R*cos(θ);
Y(R,θ) = R*sin(θ);

# Triangluar Boundary
Xₜ(R,T) = ((0<=T) & (T<=1)) * (R/2-T*R) + 
    ((1<T) & (T<=2)) * (-R/2 + (T-1)*R/2) + 
    ((2<T) & (T<=3)) * ((T-2)*R/2) 

Yₜ(R,T) = ((0<=T) & (T<=1)) * ((R*sin(π/3))/2) + 
    ((1<T) & (T<=2)) * ((R*sin(π/3))/2 - (T-1)*(R*sin(π/3))) + 
    ((2<T) & (T<=3)) * (-(R*sin(π/3))/2 + (T-2)*(R*sin(π/3)))

# Square Boundary
Xₛ(R,T) = ((0<=T) & (T<=1)) * (R*T-R/2) + 
    ((1<T) & (T<=2)) * (R/2) + 
    ((2<T) & (T<=3)) * (R/2 - R*(T-2)) + 
    ((3<T) & (T<=4)) * (-R/2);
Yₛ(R,T) = ((0<=T) & (T<=1)) * (-R/2) + 
    ((1<T) & (T<=2)) * (-R/2 + R*(T-1)) + 
    ((2<T) & (T<=3)) * (R/2) + 
    ((3<T) & (T<=4)) * (R/2 - R*(T-3));;

# Hexagon Boundary
Vertex(R) = [R R*cos(pi/3) R*cos(2*pi/3) R*cos(pi) R*cos(4*pi/3) R*cos(5*pi/3) R*cos(2*pi); 0 R*sin(pi/3) R*sin(2*pi/3) R*sin(pi) R*sin(4*pi/3) R*sin(5*pi/3) R*sin(2*pi)];

Xₕ(R,T) = ((0<=T) & (T<=1)) * (Vertex(R)[1,1] + 
    T *(Vertex(R)[1,2] - Vertex(R)[1,1])) + 
    ((1<T) & (T<=2)) * (Vertex(R)[1,2] + 
    (T-1)*(Vertex(R)[1,3] - Vertex(R)[1,2]))+ 
    ((2<T) & (T<=3)) * (Vertex(R)[1,3] + 
    (T-2)*(Vertex(R)[1,4] - Vertex(R)[1,3]))+ 
    ((3<T) & (T<=4)) * (Vertex(R)[1,4] + 
    (T-3)*(Vertex(R)[1,5] - Vertex(R)[1,4]))+ 
    ((4<T) & (T<=5)) * (Vertex(R)[1,5] + 
    (T-4)*(Vertex(R)[1,6] - Vertex(R)[1,5]))+ 
    ((5<T) & (T<6))  * (Vertex(R)[1,6] + 
    (T-5)*(Vertex(R)[1,7] - Vertex(R)[1,6]))

Yₕ(R,T) = ((0<=T) & (T<=1)) * (Vertex(R)[2,1] + 
    T*(Vertex(R)[2,2] - Vertex(R)[2,1]))+ 
    ((1<T) & (T<=2)) * (Vertex(R)[2,2] + 
    (T-1)*(Vertex(R)[2,3] - Vertex(R)[2,2]))+ 
    ((2<T) & (T<=3)) * (Vertex(R)[2,3] + 
    (T-2)*(Vertex(R)[2,4] - Vertex(R)[2,3]))+ 
    ((3<T) & (T<=4)) * (Vertex(R)[2,4] + 
    (T-3)*(Vertex(R)[2,5] - Vertex(R)[2,4]))+ 
    ((4<T) & (T<=5)) * (Vertex(R)[2,5] + 
    (T-4)*(Vertex(R)[2,6] - Vertex(R)[2,5]))+ 
    ((5<T) & (T<6))  * (Vertex(R)[2,6] + 
    (T-5)*(Vertex(R)[2,7] - Vertex(R)[2,6]))

Xᵩ(T) = T
Yᵩ(T) = 2 + 0.5*cos(3*T)


function initial_pos_1D(u0,N,η,kf,l₀)
    kₛ = 0.5*N
    η = η/N
    kf = kf/N
    Tmax = 10;
    δt = 0.01
    p = (N,kₛ,η,kf,l₀,δt)
    tspan = (0.0,Tmax)
    prob = ODEProblem(ODE_fnc_1D_init!,u0,tspan,p)
    savetimes = LinRange(0, Tmax, 2)
    init_pos = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
    return init_pos.u[2]
end

function initial_pos_2D(u0,N,η,kf,l₀)
    kₛ = 0.5*N
    η = η/N
    kf = kf/N
    Tmax = 10;
    δt = 0.0005
    p = (N,kₛ,η,kf,l₀,δt)
    tspan = (0.0,Tmax)
    prob = ODEProblem(ODE_fnc_2D_init!,u0,tspan,p)
    savetimes = LinRange(0, Tmax, 2)
    init_pos = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
    return init_pos.u[2]
end

function u0SetUp(btype,R₀,N)
    # setting up initial conditions
    θ = collect(LinRange(0.0, 2*π, N+1))  # just use collect(θ) to convert into a vector
    pop!(θ)
    u0 = ElasticArray{Float64}(undef,2,N)
    for i in 1:N
        if btype == "circle"
            R = R₀ # to produce identical areas
            @views u0[:,i] .= [X(R,θ[i]), Y(R,θ[i])];
        elseif btype == "triangle"
            R = √((2*π*R₀^2)/sin(π/3))
            @views u0[:,i] .= [Xₜ(R,θ[i]*3/(2*π)), Yₜ(R,θ[i]*3/(2*π))]
        elseif btype == "square"
            R = √(π*(R₀^2)) # to produce identical areas
            @views u0[:,i] .= [Xₛ(R,θ[i]*2/pi), Yₛ(R,θ[i]*2/pi)]
        elseif btype == "hex"
            R = √((2/3√3)*π*(R₀^2)) # to produce identical areas
            @views u0[:,i] .= [Xₕ(R,θ[i]*3/pi), Yₕ(R,θ[i]*3/pi)]
        elseif btype == "SineWave"
            θ = collect(LinRange(0.0, 2*π, N+1))  # just use collect(θ) to convert into a vector
            pop!(θ)
            @views u0[:,i] .= [Xᵩ(θ[i]), Yᵩ(θ[i])];
        end
    end

    """
    if btype == "SineWave"
        relax_pos = initial_pos_1D(u0,N,1,0,1e-3)
    else
        relax_pos = initial_pos_2D(u0,N,1,0,1e-3)
    end
    """
    relax_pos = u0
    return oftype(relax_pos, relax_pos')
end

