
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
    kₛ = 2*N
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
    kₛ = 1
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
            #pop!(θ)
            @views u0[:,i] .= [Xᵩ(θ[i]), Yᵩ(θ[i])];
        elseif btype == "star"
            star_points = 5
            r₀ = R₀/2 
            Rotation_Angle = pi/2
            rotation_angle = Rotation_Angle + pi/star_points
            verts = StarVerticies(star_points, R₀, Rotation_Angle, r₀, rotation_angle)
            u0 = interpolate_vertices(verts, Int( round(N/(2*star_points))))'
        end
    end


    if btype == "SineWave"
        relax_pos = initial_pos_1D(u0',N,1,0,1e-3)
    else
        relax_pos = initial_pos_2D(u0',N,1,0,1e-3)
    end
    return relax_pos
    """
    relax_pos = u0

    return oftype(ElasticArray{Float64}(undef,2,size(relax_pos,2)), hcat(relax_pos)')
    """
end

##############################################################################################

# linear interpolation for custom shapes

function interpolate_vertices(vertices, N)
    # Ensure at least two vertices are given
    if length(vertices) < 2
        error("At least two vertices are needed.")
    end

    interpolated_x = Float64[]
    interpolated_y = Float64[]

    for i in 1:size(vertices,1)-1
        x1, y1 = vertices[i,:]
        x2, y2 = vertices[i+1,:]

        # Calculate step sizes for x and y coordinates
        step_x = (x2 - x1) / N
        step_y = (y2 - y1) / N

        for j in 0:N-1
            x = x1 + j * step_x
            y = y1 + j * step_y
            push!(interpolated_x, x)
            push!(interpolated_y, y)
        end
    end

    return hcat([interpolated_x, interpolated_y]...)
end


function regular_polygon_vertices(N, R, rotation_angle)
    vertices = Vector{Float64}[]

    for i in 0:N-1
        angle = 2π * i / N + rotation_angle
        x = R * cos(angle)
        y = R * sin(angle)
        push!(vertices, [x, y])
    end

    return hcat(vertices...)'
end

function StarVerticies(N, R, Rotation_Angle, r, rotation_angle)
    # finding R such that areas will match with formula A = 2N(0.5*R₀*rₒ*sin(θ))
    R₀ = √((4*π*(R^2))/(2*N*sin(π/N)))
    rₒ = R₀/2
    # empty vector
    StarVerts = Vector{Float64}[]
    # generating verticies for outside polygon
    VERTS = regular_polygon_vertices(N, R₀, Rotation_Angle)
    # generating verticies for inside polygon
    verts = regular_polygon_vertices(N, rₒ, rotation_angle)

    # combining the two
    for i in 1:N
        push!(StarVerts, VERTS[i,:])
        push!(StarVerts, verts[i,:])
    end
    push!(StarVerts, VERTS[1,:])

    return hcat(StarVerts...)'
end


# Sample code for usage
"""
N = 7
R = 1
r = 0.5

Rotation_Angle = pi/2
rotation_angle = Rotation_Angle + pi/N

verts = StarVerticies(N, R, Rotation_Angle, r, rotation_angle)
scatter(verts[:,1],verts[:,2])
u = interpolate_vertices(verts, k)

scatter(u[:,1],u[:,2])
"""
