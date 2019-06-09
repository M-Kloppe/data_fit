%Funktion die zu Freiheitsgraden aus mdsps-Funktion die zugehörigen Punkte
%ermittelt

%M. Kloppe, Juni 2019

%Parameter:
%dof    ...Vektor der FG
%nv     ...Anzahl der Eckpunkte in ursprünglicher TRiangulierung
%x,y    ...Koordinaten der Punkte in PS-Verfeinerung
%nt     ...Anzahl der Dreiecke der ursprünglichen Triangulierung
%ie1o,ie2o ...Indizes der Eckpunkte, die eine Kante <ie1,ie2> bilden (urspr.
%Triang)

function [xdof,ydof]=finddofs(dof,nv,nt,x,y,ie1o,ie2o)
%Speicherreservierung
xdof=zeros(length(dof),1);
ydof=xdof;

%Anzahl der Kanten
ne=length(ie1o);

%Anzahl der möglichen Punkte in MDS
maxd=nv+nt+ne*3;

%Nummer des kleinstmögl. FG in MDS, der kein Eckpunkt ist
mind=nv+nt+ne+1;

%Hilfsvektor, der anzeigt ob entsprechender Freiheitsgrad auftritt
%1-> tritt auf, 0 sonst
helpd=zeros(maxd,1);


for i=1:length(dof)
    helpd(dof(i))=1;
end


%die ersten nv Freiheitsgrade sind die Eckpunkte der ursprünglichen
%Triangulierung

xdof(1:nv)=x(1:nv);
ydof(1:nv)=y(1:nv);

%zähler
zp=nv+1;

%schleife ermittelt nun die Punkte zu den FG
%Voraussetzung: keine Lücken in Triangulierung
for i=mind:maxd   
    if helpd(i)==1
        %ermittle Kante, auf der FG liegt
        if mod(i,2)==1
            
            kante=1/2*(i-nv-nt-ne);
            
            xdof(zp)=1/2*(x(ie2o(kante))+x(nv+nt+kante));
            ydof(zp)=1/2*(y(ie2o(kante))+y(nv+nt+kante));
            
            zp=zp+1;
        else
            kante=1/2*(i+1-nv-nt-ne);
            
            xdof(zp)=1/2*(x(ie1o(kante))+x(nv+nt+kante));
            ydof(zp)=1/2*(y(ie1o(kante))+y(nv+nt+kante));
            zp=zp+1;
        end
        
    end
end

if length(dof)~= length(xdof)
    error('Anzahl ermittelter Punkte stimmt nicht!');
end

end