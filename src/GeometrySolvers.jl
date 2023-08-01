function lineIntersection(rₘ₁,rₗ,rₘ₂,rᵣ)
    @views r = rₗ-rₘ₁
    @views s = rᵣ-rₘ₂
    
    @views d = r[1]*s[2] - r[2]*s[1]
    @views u = ((rₘ₂[1] - rₘ₁[1])*r[2] - (rₘ₂[2] - rₘ₁[2])*r[1])/d
    @views t = ((rₘ₂[1] - rₘ₁[1])*s[2] - (rₘ₂[2] - rₘ₁[2])*s[1])/d

    if(0<=u && u<=1 && 0<=t && t<=1)
        #println("Yes these intersect at: ")
        #println(rₘ₁ + t*r)
        return(rₘ₁ + t*r)
    else
        #println("No these lines dont intersect, midpoint: " )
        #println(((rₘ₁ + rₘ₂)/2))
        return (rₘ₁ + rₘ₂)/2
    end
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
end