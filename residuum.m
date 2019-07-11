%Funktion berechnet euklid. Norm und Maximumsnorm des Residuenvektors der
%Splinefunktion

%M. Kloppe, Juni 2019
function [res,maxres,r]=residuum(TRI,x,y,xd,yd,zd,c,v1,v2,v3,e1,e2,e3)
%TRI - ConnectivityList der Triangulierung
%x,y - kart. Koord. der Eckpunkte
%xd,yd,zd - kart. Koordinaten der Datenpunkte
%c - Vektor der B-Koeffizienten des Splines
% v1,v2,v3 - Liste der Eckenindizes der Dreiecke
% e1,e2,e3 - Liste der Kantenindizes der Dreiecke

%addpath('splinepak');

%erzeuge Triangulierung
T=triangulation(TRI,[x,y]);

%Initialisierung des Residuenvektors
r=zeros(length(zd),1);

%Anzahl der Eckpunkte
nv=length(x);

%Ordne alle Punkte den Dreiecken zu
ti=pointLocation(T,xd,yd);

%Ermittle baryzentrische Koordinaten der Datenpunkte
 B = cartesianToBarycentric(T,ti,[xd,yd]);

%Schleife über alle Datenpunkte
for i=1:length(zd)
    %Reduziere Spline s auf das Dreieck, in dem aktueller Punkt liegt
    %wähle dazu die zugehörigen B-Koeffizienten aus
    index=finde_ind(ti(i),nv,v1,v2,v3,e1,e2,e3);
    ci=c(index);

    %Ermittle baryzentrische Koordinaten der Datenpunkte
    %B = cartesianToBarycentric(T,ti(i),[xd(i),yd(i)]);
 
    %Werte Spline mit Hilfe des De Casteljau Algorithmus aus (verwende hier
    % Funktion aus Schumaker SplinePAK 2014)
    val=decast(2,B(i,1),B(i,2),B(i,3),ci);

    %Aktualisiere Residuenvektor
    r(i)=zd(i)-val;
end
%Maximumsnorm des Residuums
maxres=norm(r,inf);

%Euklidische Norm des Residuums
res=norm(r);




