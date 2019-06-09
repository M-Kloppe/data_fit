%Funktion berechnet Matrix MT für Assemblierungsalgorithmus der
%Energiematrix M

%M. Kloppe, Juni 2019

%x,y sind Spaltenvektoren mit den Eckpunkten des aktuellen Dreiecks
function [MT]=build_MT(x,y)
%Speicherreservierung
MT=zeros(6,6);
    
%Berechne Fläche des Dreiecks
area=polyarea(x,y);

%Berechne Transformationsmatrix T des Einheitsdreiecks
 T=[([x(2),y(2)]-[x(1),y(1)])',([x(3),y(3)]-[x(1),y(1)])'];
 
%Determinante zu T
 d=det(T);

        
%Partielle Abl. der Koordianten xi(x,y), eta(x,y) des Referenzdreiecks
 xi_x=T(2,2)/d;
 xi_y=-T(1,2)/d;
 eta_x=-T(2,1)/d;
 eta_y=T(1,1)/d;
        
%Zweite partielle Ableitungen der Bernsteinpolynome (d=2) auf
%aktuellem Dreieck als Vektoren
%Bxx ist ein Vektor, der die sechs Ableitungen DxxBi in lexikographischer 
%Reihenfolge enthält
%Analog Bxy, Byy
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