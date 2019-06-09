%Funktion zur Berechnung der Matrix H = G + \mu * M
%Funktion erstellt die Grammatrix G, die Energiematrix M sowie die rechte
%Seite d, die für die Berechnung des (gegl.) LSQ-Splines benötigt werden

%M. Kloppe, Juni 2019

% Parameter
% A - Transforationsmatrix
% x,y - Kartesische Koordinaten der Eckpunkte
% v1,v2,v3 - Liste der Eckenindizes der Dreiecke
% e1,e2,e3 - Liste der Kantenindizes der Dreiecke
%xd,yd,zd  - Datenpunkte
%wd - Gewichte
%nv - Anzahl der Eckpunkte

function [G,d,M,H]=build_H_d(A,x,y,v1,v2,v3,e1,e2,e3,xd,yd,zd,wd,lambda)
%Wir wollen den Algorithmus 7.3 aus Schumaker - Computational Methods
%verwenden

%Dimension von G = Spaltenanzahl der Trafomatrix A = Dimension des
%Splineraums
N=size(A,2);


%Initialisierung von G, M und d
G=zeros(N,N);
M=zeros(N,N);
d=zeros(N,1);


%Aktuelle Triangulierung
Tlist=[v1,v2,v3];
P=[x,y];
TRI=triangulation(Tlist,P);


%Anzahl der Dreiecke
nt=length(v1);

%Anzahl aller Eckpunkte
nv=length(x);

%Ordne alle Punkte den Dreiecken zu
ti=pointLocation(TRI,xd,yd);


%Assemblierung von G, M und r (Schleife über alle Dreiecke)
for i=1:nt
    %Ermittle alle Punkte, die sich in betrachtetem Dreieck befinden
    ai=find(ti==i);

    if isempty(ai)==0
    %Erstelle die Matrix G(T), sowie r(T)
    [GT,dT]=build_GT(TRI,i,wd(ai),xd(ai),yd(ai),zd(ai));

    %Ermittle aktuelle Knotenindizes
    [index]=finde_ind(i,nv,v1,v2,v3,e1,e2,e3);

    %Wähle Untermatrix A(T) aus A aus (Zeilen zu index)
    AT=A(index,:);

    %Aktualisiere G
    G=G+AT'*GT*AT;

    %Aktualisiere r
    d=d+AT'*dT;
    
    %Ermittle MT und aktualisiere M
    if lambda>0
        vx=[x(v1(i)),x(v2(i)),x(v3(i))]';
        vy=[y(v1(i)),y(v2(i)),y(v3(i))]';
        MT=build_MT(vx,vy);
        M=M+AT'*MT*AT;
    end
    end
end
H=G+lambda*M;
end