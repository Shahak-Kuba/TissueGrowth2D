"""
function findIntersection(rₗ,rₘ₁,rₘ₂,rᵣ)
   Lₗ = Line(tuple(rₘ₁...), tuple(rₗ...))
   Lᵣ = Line(tuple(rₘ₂...), tuple(rᵣ...))
   # find intersection point
   i = tuple(coordinates(intersect(Lₗ, Lᵣ))...)
   # checking if intersection is within 

end
"""

function lineIntersection(rₘ₁,rₗ,rₘ₂,rᵣ)
    r = rₗ-rₘ₁
    s = rᵣ-rₘ₂
    
    d = r[1]*s[2] - r[2]*s[1]
    u = ((rₘ₂[1] - rₘ₁[1])*r[2] - (rₘ₂[2] - rₘ₁[2])*r[1])/d
    t = ((rₘ₂[1] - rₘ₁[1])*s[2] - (rₘ₂[2] - rₘ₁[2])*s[1])/d

    if(0<=u && u<=1 && 0<=t && t<=1)
        println("Yes these intersect at: ")
        println(rₘ₁ + t*r)
    else
        println("No these lines dont intersect, midpoint: " )
        println(((rₘ₁ + rₘ₂)/2))
    end
end

