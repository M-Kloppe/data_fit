%function determines k_th order star of one or more vertices

%M. Kloppe, Juni 2019

%Paramter:
%vlist - cartesian coordinates of the vertices
%k - order of star
%indv - Indizes of vertices
%tri - ConnectivityList of triangulation

%return values:
%vstar - coordinates of vertices in star
%indstar - Indizes of these vertices
%tristar - triangles in star
%indtri - indizes of these triangles (regarding used triangulation)
function [vstar,indstar,tristar,indtri]=starvk(vlist,k,indv,tri)
%determine star recursively
    if k==1
        vstar=[];
        indstar=[];
        tristar=[];
        indtri=[];
        for j=1:length(indv)
            %determine for each vertex all triangles containing this vertex
            m1=find(tri(:,1)==indv(j));
            m2=find(tri(:,2)==indv(j));
            m3=find(tri(:,3)==indv(j));
        
            %just use every triangle once -> all triangles in first order
            %star
            indhelp=unique([m1;m2;m3]);
            indhelptri=indhelp;
            indhelp=tri(indhelp,:);
            %just use every vertex once -> all vertices in star
            indhelp2=indhelp;
            indhelp=unique(indhelp);
        
            %update lists of vertices and triangles(+ Indizes)
            vstar=[vstar;vlist(indhelp,:)];
            indstar=[indstar;indhelp];
            tristar=[tristar;indhelp2];
            indtri=[indtri;indhelptri];
        end
    else
        indstar=[];
        tristar=[];
        indtri=[];
        %higher order star -> build star of the boundary nodes of the star
        % order one less
        for i=1:k
            [vstarhelp,indstarhelp,tristarhelp,indhelptri]=starvk(vlist,1,indv,tri);
            
            %determine boundary nodes
            bound=boundary(vstarhelp(:,1),vstarhelp(:,2),0);
            indv=indstarhelp(bound);
            
            %update Lists
            indstar=[indstar;indstarhelp];
            tristar=[tristar;tristarhelp];
            indtri=[indtri;indhelptri];
        end
        %just use vertices and triangles once
        indstar=unique(indstar);
        vstar=vlist(indstar,:);
        tristar=unique(tristar,'rows');
        indtri=unique(indtri);
    end
end