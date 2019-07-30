%function computes matrix MT for assembling algorithm of energy matrix M

%M. Kloppe, Juni 2019

%x,y column vectors with vertices of current triangle
function [MT]=build_MT(x,y)
%Initialise
MT=zeros(6,6);
    
%compute surface of triangle
area=polyarea(x,y);

%comute transformationmatrix T of unit triangle
 T=[([x(2),y(2)]-[x(1),y(1)])',([x(3),y(3)]-[x(1),y(1)])'];
 
%Determinant of T
 d=det(T);

        
%partiell derivative of coordinats xi(x,y), eta(x,y) of reference triangle
 xi_x=T(2,2)/d;
 xi_y=-T(1,2)/d;
 eta_x=-T(2,1)/d;
 eta_y=T(1,1)/d;
        
%second partiell derivative of Bernsteinpolynomials (d=2) on current
%triangle (as vectors)
%Bxx is a vector with all the six  derivatives DxxBi in lexikographical order
%analogue Bxy, Byy
Bxx=[2*xi_x^2;...
    4*xi_x*eta_x;...
    -4*xi_x^2-4*xi_x*eta_x;...
    2*eta_x^2;...
    -4*xi_x*eta_x-4*eta_x^2;...
    2*xi_x^2+4*xi_x*eta_x+2*eta_x^2];

Bxy=[2*xi_x*xi_y;...
    2*xi_x*eta_y+2*xi_y*eta_x;...
    -4*xi_x*xi_y-2*xi_x*eta_y-2*xi_y*eta_x;...
    2*eta_x*eta_y;...
    -2*xi_x*eta_y-2*xi_y*eta_x-4*eta_x*eta_y;...
    2*xi_x*xi_y+2*xi_x*eta_y+2*xi_y*eta_x+2*eta_x*eta_y];

Byy=[2*xi_y^2;...
    4*xi_y*eta_y;...
    -4*xi_y^2-4*xi_y*eta_y;...
    2*eta_y^2;...
    -4*xi_y*eta_y-4*eta_y^2;...
    2*xi_y^2+4*xi_y*eta_y+2*eta_y^2];


        
%Matrix MT
 MT=area*((Bxx*Bxx')+2*(Bxy*Bxy')+(Byy*Byy'));
end