%Funktion wandelt gegebenen Höhendaten (Matrix B - Stromboli) in passendes
%Format um, d.h. Fkt. erzeugt Vektoren xd,yd,zd, die kartesische 
%Koordinaten der Datenpunkte enthalten

%M. Kloppe, Juni 2019
function [xd,yd,zd,nd]=create_data(B)
%Dimension der Matrix B
[z,s]=size(B);

%Hilfsvektoren zur Erzeugung von xd und yd
xhelp=linspace(1,s,s);
yhelp=linspace(1,z,z);

%hier entsprechen xd und yd "Pixeln" in der xy-Ebene
%Erzeuge uniformes Gitter und Ordne die Elemente als Vektoren an
[X,Y]=meshgrid(xhelp,yhelp);
xd=reshape(X,[length(xhelp)*length(yhelp),1]);
yd=reshape(Y,[length(xhelp)*length(yhelp),1]);

%Erzeuge Vektor der Höhendaten
zd=reshape(B,[length(xhelp)*length(yhelp),1]);
zd=double(zd);
nd=length(zd);
end