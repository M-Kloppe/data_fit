%Funktion gibt für Dreieck mit den Eckpunkten x,y,z (x = linke untere Ecke) 
%die Knoten für die quadratische Interpolation mit Lagrange-El. zurück

% außerdem werden Transformationsmatrix T und die Verschiebung b bzgl des Eiheitsdreiecks ausgegeben 

function [x1,x2,x3,x4,x5,x6,T,b]=quad_nodes_lexi(x,y,z)
    x1=x;
    x6=z;
    x4=y;
    
    x3=[1/2*(x(1)+z(1)),1/2*(x(2)+z(2))];
    x2=[1/2*(x(1)+y(1)),1/2*(x(2)+y(2))];
    x5=[1/2*(z(1)+y(1)),1/2*(z(2)+y(2))];
    
    T=[(y-x)',(z-x)'];
    b=x;
end