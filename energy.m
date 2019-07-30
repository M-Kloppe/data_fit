%Function computes the energy of a PS-Spline with coefficient vector c

%M. Kloppe, Juni 2019

% x,y - coordinates of the vertices of the triangulation
function [E]=energy(x,y,c,v1,v2,v3,e1,e2,e3)
    
%coordinates of the vertices of each triangle
vx=[x(v1),x(v2),x(v3)]';
vy=[y(v1),y(v2),y(v3)]';
    
%area of all the triangles
area=polyarea(vx,vy);

%number of triangles
nt=length(v1);
    
%number of vertices
nv=length(x);
    
%summing all the energies on the triangles in E up
E=0;
    
%loop over all triangles
for i=1:nt
 %transformation matrix
  T=[([x(v2(i)),y(v2(i))]-[x(v1(i)),y(v1(i))])',...
      ([x(v3(i)),y(v3(i))]-[x(v1(i)),y(v1(i))])'];
     
  %Determinant of T
  D=det(T);

  %partiell derivatives of the coordinates xi(x,y) eta(x,y) of the
  %reference triangle
  xi_x=T(2,2)/D;
  xi_y=-T(1,2)/D;
  eta_x=-T(2,1)/D;
  eta_y=T(1,1)/D;

%second partiell derivative of Bernsteinpolynomials (d=2) on current
%triangle (as vectors)
%Bxx is a vector with all the six  derivatives DxxBi in lexikographical order
%analogue Bxy, Byy
   Bxx=[2*xi_x^2;2*eta_x^2;2*xi_x^2+4*xi_x*eta_x+2*eta_x^2;4*xi_x*eta_x;...
       -4*xi_x^2-4*xi_x*eta_x;-4*xi_x*eta_x-4*eta_x^2];
   Bxy=[2*xi_x*xi_y;2*eta_x*eta_y; 2*xi_x*xi_y+2*xi_x*eta_y+...
       2*xi_y*eta_x+2*eta_x*eta_y;2*xi_x*eta_y+2*xi_y*eta_x;...
       -4*xi_x*xi_y-2*xi_x*eta_y-2*xi_y*eta_x;-2*xi_x*eta_y-...
       2*xi_y*eta_x-4*eta_x*eta_y];
   Byy=[2*xi_y^2;2*eta_y^2;2*xi_y^2+4*xi_y*eta_y+2*eta_y^2;4*xi_y*eta_y;...
       -4*xi_y^2-4*xi_y*eta_y;-4*xi_y*eta_y-4*eta_y^2];
        
  %find the indizes of the nodes on the current triangle 
  [index]=finde_ind(i,nv,v1,v2,v3,e1,e2,e3);
  ci=c(index);

  %energy of the spline on the current triangle
   Ei=area(i)*((ci'*Bxx)^2+2*(ci'*Bxy)^2+(ci'*Byy)^2);

  %update E
   E=E+Ei;
end
end