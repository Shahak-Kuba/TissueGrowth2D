function lineIntersection(rₘ₁,rₗ,rₘ₂,rᵣ)
    intersect = zeros(size(rₘ₁))
    r = rₗ.-rₘ₁
    s = rᵣ.-rₘ₂
    
    d = r[:,1].*s[:,2] - r[:,2].*s[:,1]
    # performing determinant test in case lines are parallel
    for i in 1:size(d,1)
        if d[i] == 0
            intersect[i,:] = (rₘ₁[i,:] + rₘ₂[i,:])./2
        
        else
            u = ((rₘ₂[i,1] - rₘ₁[i,1])*r[i,2] - (rₘ₂[i,2] - rₘ₁[i,2])*r[i,1])/d[i]
            t = ((rₘ₂[i,1] - rₘ₁[i,1])*s[i,2] - (rₘ₂[i,2] - rₘ₁[i,2])*s[i,1])/d[i]

            if 0≤u≤1 && 0≤t≤1
                #println("Yes these intersect at: ")
                #println(rₘ₁ + t*r)
                intersect[i,:] = (rₘ₁[i,:] + t*r[i,:])
            else
                #println("No these lines dont intersect, midpoint: " )
                #println(((rₘ₁ + rₘ₂)/2))
                intersect[i,:] =  (rₘ₁[i,:] + rₘ₂[i,:])/2
            end
        end
    end
    return intersect
end

function Ω(p)
    A = 0
    for ii in axes(p,2)
        if ii == size(p,2)
            A += (p[ii,1]*p[1,2] -  p[ii,2]*p[1,1])
        else
            A += (p[ii,1]*p[ii+1,2] -  p[ii,2]*p[ii+1,1])
        end
    end
    return abs(A)/2;
end


using Test
# running test cases
@testset "LineIntersections" begin
    p1 = [0;0]
    p2 = [0;1]
    p3 = [1;0]
    p4 = [√2;√2]

    @test lineIntersection(p1,p4,p2,p3) ≈ [0.5;0.5] atol=0.01
    @test lineIntersection(p1,p3,p2,p4) ≈ [0.;0.5] atol=0.01
    @test lineIntersection(p1,p2,p3,p4) ≈ [0.5;0.] atol=0.01
    @test lineIntersection(p1,p2,p1,p3) ≈ [0.;0.] atol=0.01
    @test lineIntersection(p3,p1,p2,p4) ≈ [0.5;0.5] atol=0.01
    @test lineIntersection(p1,p1,p2,p3) ≈ [0.;0.5] atol=0.01    # det(M) == 0 test
end

#@testset "ShapeArea" begin
    # square geometry
    p = [0 0 1 1; 0 1 1 0]
    p2 = [0 0 2 2; 0 2 2 0]

    @test Ω(p) ≈ 1 atol=0.01
    @test Ω(p2) ≈ 4 atol=0.01

    N = 288
    R₀ = 2

    #testing circle
    btype = "circle"
    P = u0SetUp(btype,R₀,N)
    @test Ω(P) ≈ 2*π*R₀ atol=0.01

    """
    btype = "hex"
    Pₕ = u0SetUp(btype,R₀,N)
    @test Ω(Pₕ) ≈ 3*√3 / 2 * (√((2/3√3)*π*(R₀^2)))^2 atol=0.01

    btype = "triangle"
    Pₜ = u0SetUp(btype,R₀,N)
    @test Ω(Pₜ) ≈ √3 / 4 * (√((2*π*R₀^2)/sin(π/3)))^2 atol=0.01
    """

#end