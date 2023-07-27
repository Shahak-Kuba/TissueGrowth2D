
# Circular Boundary
X(R,θ) = R*cos(θ);
Y(R,θ) = R*sin(θ);

# Square Boundary
#Xₛ(R,T) = ((0<=T) & (T<=1)) * (R) + ((1<T) & (T<=2)) * (R - 0.1 + (0.1 * cos(0.5*pi*(T-1)))) + ((2<T) & (T<=3)) * (-(R - 0.1)*(2*T-5)) + ((3<T) & (T<=4)) * (-R + 0.1 - (0.1 * sin(0.5*pi*(T-3)))) + ((4<T) & (T<=5)) * (-R) + ((5<T) & (T<=6)) * (-R + 0.1 - (0.1 * cos(0.5*pi*(T-5)))) + ((6<T) & (T<=7)) * ((R-0.1)*(2*T - 13)) + ((7<T) & (T<=8)) * (R - 0.1 + (0.1*sin(0.5*pi*(T-7))));
#Yₛ(R,T) = ((0<=T) & (T<=1)) * (-(R - 0.1)*(2*T-1)) + ((1<T) & (T<=2)) * (-R + 0.1 - 0.1 * sin(0.5*pi*(T-1))) + ((2<T) & (T<=3)) * (-R) + ((3<T) & (T<=4)) * (-R + 0.1 - 0.1 * cos(0.5*pi*(T-3))) + ((4<T) & (T<=5)) * ((R - 0.1)*(2*T - 9)) + ((5<T) & (T<=6)) * (R - 0.1 + 0.1 * sin(0.5*pi*(T-5))) + ((6<T) & (T<=7)) * (R) + ((7<T) & (T<=8)) * (R - 0.1 + 0.1 * cos(0.5*pi*(T-7)));

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
Xₕ(R,T) = ((0<=T) & (T<=1)) * (Vertex(R)[1,1] + T *(Vertex(R)[1,2] - Vertex(R)[1,1])) + ((1<T) & (T<=2)) * (Vertex(R)[1,2] + (T-1)*(Vertex(R)[1,3] - Vertex(R)[1,2]))+ ((2<T) & (T<=3)) * (Vertex(R)[1,3] + (T-2)*(Vertex(R)[1,4] - Vertex(R)[1,3]))+ ((3<T) & (T<=4)) * (Vertex(R)[1,4] + (T-3)*(Vertex(R)[1,5] - Vertex(R)[1,4]))+ ((4<T) & (T<=5)) * (Vertex(R)[1,5] + (T-4)*(Vertex(R)[1,6] - Vertex(R)[1,5]))+ ((5<T) & (T<6))  * (Vertex(R)[1,6] + (T-5)*(Vertex(R)[1,7] - Vertex(R)[1,6]))
Yₕ(R,T) = ((0<=T) & (T<=1)) * (Vertex(R)[2,1] + T*(Vertex(R)[2,2] - Vertex(R)[2,1]))+ ((1<T) & (T<=2)) * (Vertex(R)[2,2] + (T-1)*(Vertex(R)[2,3] - Vertex(R)[2,2]))+ ((2<T) & (T<=3)) * (Vertex(R)[2,3] + (T-2)*(Vertex(R)[2,4] - Vertex(R)[2,3]))+ ((3<T) & (T<=4)) * (Vertex(R)[2,4] + (T-3)*(Vertex(R)[2,5] - Vertex(R)[2,4]))+ ((4<T) & (T<=5)) * (Vertex(R)[2,5] + (T-4)*(Vertex(R)[2,6] - Vertex(R)[2,5]))+ ((5<T) & (T<6))  * (Vertex(R)[2,6] + (T-5)*(Vertex(R)[2,7] - Vertex(R)[2,6]))