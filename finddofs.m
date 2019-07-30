%function determine for degress of freedom related cartesian coordinates

%M. Kloppe, Juni 2019

%Parameter:
%dof    ...vector of degrees of freedom
%nv     ...number of vertices in triangulation TRIo
%x,y    ...coordinates of points in PS-refinement
%nt     ...number of triangles in triangulation
%ie1o,ie2o ...Indizes of vertices on one edge <ie1,ie2> (in TRIo)

function [xdof,ydof]=finddofs(dof,nv,nt,x,y,ie1o,ie2o)
%initialize
xdof=zeros(length(dof),1);
ydof=xdof;

%number of edges
ne=length(ie1o);

%(maximal) number of points in MDS
maxd=nv+nt+ne*3;

%index of the smallest possible dof in MDS (no vertex)
mind=nv+nt+ne+1;

%flag for degrees of freedom
%1-> dof is used, 0 otherwise
helpd=zeros(maxd,1);


for i=1:length(dof)
    helpd(dof(i))=1;
end


%first nv dofs are the vertices of TRIo

xdof(1:nv)=x(1:nv);
ydof(1:nv)=y(1:nv);

%counter
zp=nv+1;

%loop determines points for the dofs
%assumption: no holes in triangulation
for i=mind:maxd   
    if helpd(i)==1
        %determine edge with one dof
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