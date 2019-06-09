%Funktion berechnet die Energie eines PS-Splines mit B-Koeffis c

%M. Kloppe, Juni 2019

% x,y - Koordinaten der Eckpunkte der Triangulierung
function [E]=energy(x,y,c,v1,v2,v3,e1,e2,e3)
    
%Erzeuge Matrizen mit den Koordinaten der Eckpunkte jedes Dreiecks
vx=[x(v1),x(v2),x(v3)]';
vy=[y(v1),y(v2),y(v3)]';
    
%Berechne Fläche aller Dreiecke
area=polyarea(vx,vy);

%Anzahl der Dreiecke
nt=length(v1);
    
%Anzahl der Eckpunkte
nv=length(x);
    
%E entsteht durch Summierung der "Energien" auf den einzelnen Dreiecken
E=0;
    
%Schleife über alle Dreiecke
for i=1:nt
 %Berechne Transformationsmatrix T des Einheitsdreiecks
  T=[([x(v2(i)),y(v2(i))]-[x(v1(i)),y(v1(i))])',...
      ([x(v3(i)),y(v3(i))]-[x(v1(i)),y(v1(i))])'];
     
  %Determinante zu T
  D=det(T);

  %Partielle Ableitungen der Koordianten xi(x,y) eta(x,y) des
  %Referenzdreiecks
  xi_x=T(2,2)/D;
  xi_y=-T(1,2)/D;
  eta_x=-T(2,1)/D;
  eta_y=T(1,1)/D;

  %Zweite partielle Ableitungen der Bernsteinpolynome (d=2) auf
  %aktuellem Dreieck als Vektoren
  %Bxx ist ein Vektor, der die sechs Ableitungen DxxBi enthält
  %Analog Bxy, Byy
   Bxx=[2*xi_x^2;2*eta_x^2;2*xi_x^2+4*xi_x*eta_x+2*eta_x^2;4*xi_x*eta_x;...
       -4*xi_x^2-4*xi_x*eta_x;-4*xi_x*eta_x-4*eta_x^2];
   Bxy=[2*xi_x*xi_y;2*eta_x*eta_y; 2*xi_x*xi_y+2*xi_x*eta_y+...
       2*xi_y*eta_x+2*eta_x*eta_y;2*xi_x*eta_y+2*xi_y*eta_x;...
       -4*xi_x*xi_y-2*xi_x*eta_y-2*xi_y*eta_x;-2*xi_x*eta_y-...
       2*xi_y*eta_x-4*eta_x*eta_y];
   Byy=[2*xi_y^2;2*eta_y^2;2*xi_y^2+4*xi_y*eta_y+2*eta_y^2;4*xi_y*eta_y;...
       -4*xi_y^2-4*xi_y*eta_y;-4*xi_y*eta_y-4*eta_y^2];
        
  %Wähle B-Koeffis zu Knoten auf Dreieck
  [index]=finde_ind(i,nv,v1,v2,v3,e1,e2,e3);
  ci=c(index);

  %Energie des Splines auf aktuellem Dreieck
   Ei=area(i)*((ci'*Bxx)^2+2*(ci'*Bxy)^2+(ci'*Byy)^2);

  %Aktualisierung der Summe der Energien
   E=E+Ei;
end
end