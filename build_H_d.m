%function computes matrix H = G + \mu * M
%therefor get Grammatrix G, Energymatrix M and right-hand-side vector d
%for solving the LSQ-Splines-equation

%M. Kloppe, Juni 2019

% Parameter
% A - Transforationmatrix
% x,y - Cartesian koordinates of the vertices
% v1,v2,v3 - list of indizes of vertices (each triangle)
% e1,e2,e3 - list of indizes of edges (each triangle)
%xd,yd,zd  - data points
%wd - weights
%nv - number of vertices

function [G,d,M,H]=build_H_d(A,x,y,v1,v2,v3,e1,e2,e3,xd,yd,zd,wd,mu)
%we are using algorithm 7.3 from Schumaker - Computational Methods

%Dimension of G = number of columnsof Trafomatrix A = Dimension of
%Splinespace
N=size(A,2);


%Initialise of G, M and d
G=zeros(N,N);
M=zeros(N,N);
d=zeros(N,1);


%current triangulation
Tlist=[v1,v2,v3];
P=[x,y];
TRI=triangulation(Tlist,P);


%number of triangles
nt=length(v1);

%number of vertices
nv=length(x);

%find for each point the triangle
ti=pointLocation(TRI,xd,yd);


%Assembly of G, M and d (loop over all triangles)
for i=1:nt
    %determine all points in current triangle
    ai=find(ti==i);

    if isempty(ai)==0
    %build G(T) and r(T)
    [GT,dT]=build_GT(TRI,i,wd(ai),xd(ai),yd(ai),zd(ai));

    %determine current indizes of the nodes
    [index]=finde_ind(i,nv,v1,v2,v3,e1,e2,e3);

    %Choose submatrix A(T) of A (rows -> index)
    AT=A(index,:);

    %Update G
    G=G+AT'*GT*AT;

    %Update d
    d=d+AT'*dT;
    
    %determine MT and update M
    if mu>0
        vx=[x(v1(i)),x(v2(i)),x(v3(i))]';
        vy=[y(v1(i)),y(v2(i)),y(v3(i))]';
        MT=build_MT(vx,vy);
        M=M+AT'*MT*AT;
    end
    end
end
H=G+mu*M;
end