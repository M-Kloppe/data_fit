%Funktion ermittelt Stern k-ter Ordnung eines (bzw. mehrerer) Eckpunktes

%M. Kloppe, Juni 2019

%Paramter:
%vlist - Kartesische Koord. der Eckpunkte
%k - Ordnung des Sterns
%indv - Indizes der betrachteten Eckpunkte
%tri - ConnectivityList der betrachteten Triangulierung

%Rückgabewerte:
%vstar - Koordinaten der Eckpunkte in star
%indstar - zugehörige Indizes
%tristar - Dreiecke in star
%indtri - indizes der Dreiecke in betrachteter Triangulierung
function [vstar,indstar,tristar,indtri]=starvk(vlist,k,indv,tri)
%ermittle Stern rekursiv
    if k==1
        vstar=[];
        indstar=[];
        tristar=[];
        indtri=[];
        for j=1:length(indv)
            %ermittle für jeden Eckpunkte alle Dreiecke mit diesem Eckpunkt
            m1=find(tri(:,1)==indv(j));
            m2=find(tri(:,2)==indv(j));
            m3=find(tri(:,3)==indv(j));
        
            %Streiche mehrfache Dreiecke -> das sind alle Dreiecke im Stern
            indhelp=unique([m1;m2;m3]);
            indhelptri=indhelp;
            indhelp=tri(indhelp,:);
            %streiche auch mehrfache Ecken -> das sind alle Eckpkte im
            %Stern
            indhelp2=indhelp;
            indhelp=unique(indhelp);
        
            %Ergänze Liste aller Eckpunkte und Dreiecke (+zugeh. Indizes)
            vstar=[vstar;vlist(indhelp,:)];
            indstar=[indstar;indhelp];
            tristar=[tristar;indhelp2];
            indtri=[indtri;indhelptri];
        end
    else
        indstar=[];
        tristar=[];
        indtri=[];
        %für Stern höherer Ordnung erzeuge Stern der Randknoten des Sterns
        %niedrigerer Ordnung
        for i=1:k
            [vstarhelp,indstarhelp,tristarhelp,indhelptri]=starvk(vlist,1,indv,tri);
            
            %ermittle neue Randknoten des Sterns
            bound=boundary(vstarhelp(:,1),vstarhelp(:,2),0);
            indv=indstarhelp(bound);
            
            %ergänze Listen
            indstar=[indstar;indstarhelp];
            tristar=[tristar;tristarhelp];
            indtri=[indtri;indhelptri];
        end
        %streiche doppelte Knoten bzw. Dreiecke
        indstar=unique(indstar);
        vstar=vlist(indstar,:);
        tristar=unique(tristar,'rows');
        indtri=unique(indtri);
    end
end