%Funktion zur Ermittlung der Beobachtungsmatrix O des
%Kleinste-Quadrate-Problems

%M. Kloppe, Juni 2019

%Parameter:
%dof - Vektor der Freiheitsgrade des Splines
%A - Transformationsmatrix
%x,y - kartesische Koord. der Eckpunkte der verfeinerten Triangulierung
%v1,v2,v3, e1,e2,e3 - Vektoren der Ecken- und Kantenindizes jedes Dreiecks
%xo,yo - kartesische Koord. der Eckpunkte der urspr. Triangulierung
%TRIolist - Connectivity List der Eckpunkte der urspr. Triangulierung
%xd,yd - Vektoren der Daten
%ie1o,ie2o - Vektoren mit den Anfangs- bzw. Endknoten jeder Kante
function [O]=build_O(dof,A,x,y,v1,v2,v3,e1,e2,e3,xo,yo,TRIolist,xd,yd,ie1o,ie2o)
%Trafomatrix enthält in den Spalten die B-Koeffis zu den M-Basis-Splines 
A=full(A);

%Dimension des Spline-Raums
dim=length(dof);

%Initialisierung von O
O=zeros(length(xd),dim);

%urspruengliche Triangulierung
TRIo=triangulation(TRIolist,[xo,yo]);

%Anzahl der Eckpunkte (urspruengliche Triangulierung)
nvo=length(xo);

%Anzahl der Dreiecke (urspruengliche Triangulierung)
nto=size(TRIolist,1);

%Anzahl der Kanten (urspruengliche Triangulierung)
neo=length(ie1o);

%Ordne alle Punkte den Dreiecken zu
tio=pointLocation(TRIo,xd,yd);

%betrachtete Triangulierung (verfeinert)
TRIlist=[v1,v2,v3];
TRI=triangulation(TRIlist,[x,y]);
%triplot(TRI);

%Anzahl Eckpunkte (neue Triangulierung)
nv=length(x);

%Ordne alle Punkte den Dreiecken zu
ti=pointLocation(TRI,xd,yd);

for i=1:dim
    indv=dof(i);
    %ermittle aktuellen Eckpunkt für star
    if i>nvo
        if mod(dof(i),2)==1
            inde=1/2*(dof(i)-nvo-nto-neo); 
            indv=ie2o(inde);
        else
            inde=1/2*(dof(i)+1-nvo-nto-neo);
            indv=ie1o(inde);
        end
        
    end

    %Basis-Spline hat support auf star(v) (bzgl. urspr. Triangulierung)
    [~,~,~,indtrio]=starvk([xo,yo],1,indv,TRIolist);
    
    %Datenpunkte in diesen Dreiecken
    indpoio=ismember(tio,indtrio);

    %Verfeinerungsdreiecke mit diesen Daten
    indT=unique(ti(indpoio));
    
    for j=1:length(indT)
           %Punkte in aktuellem Dreieck
           indpoi=find(ti==indT(j));
           poi=[xd(indpoi),yd(indpoi)];
           
           %baryzentrische Koordinaten dieser Punkte
           B = cartesianToBarycentric(TRI,indT(j)*ones(length(indpoi),1),poi);
           
           %B-Koeffis des Basis-Splines
           index=finde_ind(indT(j),nv,v1,v2,v3,e1,e2,e3);
           cakt=A(index,i);           
           
           %ergänze O
           O(indpoi,i)=DeCasteljau(B(:,1),B(:,2),B(:,3),cakt');
    end     
end
end
