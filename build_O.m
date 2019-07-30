%function determines observationmatrix O of LSQ-Problem

%M. Kloppe, Juni 2019

%Parameter:
%dof - Vector of dofs
%A - transformationmatrix
%x,y - cart. Koord. of vertices of refined triangulation
%v1,v2,v3, e1,e2,e3 - Vectors of indizes of vertices and edges (each
%triangle)
%xo,yo - cart. Koord. of vertices of TRIo
%TRIolist - Connectivity List of TRIo
%xd,yd - vector of data
%ie1o,ie2o - Vector of start and end node (each edge)
function [O]=build_O(dof,A,x,y,v1,v2,v3,e1,e2,e3,xo,yo,TRIolist,xd,yd,ie1o,ie2o)
%Trafomatrix contains in each column B-coeffis of M-Basis-Splines 
A=full(A);

%Dimension of Spline-Space
dim=length(dof);

%Initialise of O
O=zeros(length(xd),dim);

%TRIo
TRIo=triangulation(TRIolist,[xo,yo]);

%number of vertices (TRIo)
nvo=length(xo);

%number of triangles (TRIo)
nto=size(TRIolist,1);

%number of edges (edges)
neo=length(ie1o);

%find for all data points the triangle in TRIo
tio=pointLocation(TRIo,xd,yd);

%refined triangulation TRI
TRIlist=[v1,v2,v3];
TRI=triangulation(TRIlist,[x,y]);
%triplot(TRI);

%number of vertices (TRI)
nv=length(x);

%find for all data points the triangle in TRI
ti=pointLocation(TRI,xd,yd);

for i=1:dim
    indv=dof(i);
    %determine current vertex for star
    if i>nvo
        if mod(dof(i),2)==1
            inde=1/2*(dof(i)-nvo-nto-neo); 
            indv=ie2o(inde);
        else
            inde=1/2*(dof(i)+1-nvo-nto-neo);
            indv=ie1o(inde);
        end
        
    end

    %Basis-Spline has support on star(v) (TRIo)
    [~,~,~,indtrio]=starvk([xo,yo],1,indv,TRIolist);
    
    %Data points in this triangle
    indpoio=ismember(tio,indtrio);

    %refinement triangles with these data points
    indT=unique(ti(indpoio));
    
    for j=1:length(indT)
           %points in current triangle
           indpoi=find(ti==indT(j));
           poi=[xd(indpoi),yd(indpoi)];
           
           %barycentric coordinates of these points
           B = cartesianToBarycentric(TRI,indT(j)*ones(length(indpoi),1),poi);
           
           %B-Coeffis of Basis-Splines
           index=finde_ind(indT(j),nv,v1,v2,v3,e1,e2,e3);
           cakt=A(index,i);           
           
           %update O
           O(indpoi,i)=DeCasteljau(B(:,1),B(:,2),B(:,3),cakt');
    end     
end
end
