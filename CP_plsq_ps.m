% Script for computing LSQ-Spline for scattered data

% M. Kloppe, Juni 2019

%Procedure: (1) Reading geographical data (e.g. Stromboli/Black Forest)
%           (2) Choosing a triangulation of the domain
%           (3) Computing triangulation information, MDS and trafomatrix for 
%               Powell-Sabin-Space (Using SplinePak of L.Schumaker)
%               -> SplinePAK: Copyright Larry Schumaker 2014
%           (4) Computing coefficient vector of LSQ-Spline
%           (5) Plot computed Spline-Function (Schumaker-Software)
%           (6) Computing some error measures


clear all;
%addpath('splinepak');

%Choose one data set for reading
%set = 1 -> Stromboli data
%set = 2 -> Black Forest data
set=2;

%Reading the data and if needed transform into the right shape
%xd,yd,zd - Cart. Coordinates of geographical data as columns
%nd -  number of data points
if set==1
    daten=load('stromboli.mat');
    [xd,yd,zd,nd]=create_data(daten.B);
elseif set==2
    daten=load('schwarzwald.mat');
    nd=daten.nd;
    xd=daten.xd;
    yd=daten.yd;
    zd=daten.zd;
end

%domain (->rectangle)
epsilon=eps;
k=4;
x1=min(xd)-epsilon*10^k;
y1=min(yd)-epsilon*10^k;
x2=max(xd)+epsilon*10^k;
y2=max(yd)+epsilon*10^k;

% choosing a precomputed triangulation (PDE-Toolbox)
% t=1 coarse triangulation (48 triangles)
% t=2 triangulation (130 triangles)
% t=3 fine triangulation (330 triangles)
t=2;

if t==0
    xo=[x1,x2,x1,x2]';
    yo=[y1,y1,y2,y2]';
    no=2;
    TRIo=[1,2,3;2,4,3];
end

if t==1
%reading first triangluation
Tdat=load('PDEtri1.mat');
TRIo=Tdat.TRIo;
xo = Tdat.xo;
yo = Tdat.yo;
no=length(xo);

%transform to the domain rectangle
xo=(xo-min(xo))/(max(xo)-min(xo))*abs(x2-x1)+x1;
yo=(yo-min(yo))/(max(yo)-min(yo))*abs(y2-y1)+y1;

elseif t==2  
%reading second triangluation
Tdat=load('PDEtri2.mat');
TRIo=Tdat.TRIo;
xo = Tdat.xo;
yo = Tdat.yo;
no=length(xo);

%transform to the domain rectangle
xo=(xo-min(xo))/(max(xo)-min(xo))*abs(x2-x1)+x1;
yo=(yo-min(yo))/(max(yo)-min(yo))*abs(y2-y1)+y1;

elseif t==3
%reading third triangluation
Tdat=load('PDEtri3.mat');
TRIo=Tdat.TRIo;
xo = Tdat.xo;
yo = Tdat.yo;
no=length(xo);

%transform to the domain rectangle
xo=(xo-min(xo))/(max(xo)-min(xo))*abs(x2-x1)+x1;
yo=(yo-min(yo))/(max(yo)-min(yo))*abs(y2-y1)+y1;
end

% Compute the triangulation information (Schumaker SplinePAK)
[nbo,neo,nto,v1o,v2o,v3o,e1o,e2o,e3o,ie1o,ie2o,trilo,triro,bdyo,...
   vadjo,eadjo,adjstarto,tadjo,tstarto,areao,TRIo] = trilists(xo,yo,TRIo);

 figure;
 triplot(TRIo,xo,yo);
 axis equal;

wd = ones(nd,1);

% simulation of noise
e = 0;
if e > 0
 %create new noise (Gauß)
 %rv=-1+2.*rand(nd,1);
 %zd = zd + e*rv;
 
 %read a precomputed noise vector
 noise=load('rauschen2.mat');
 zd = zd + e*noise.rv;
end

%Find a minimal determing set and transformation matrix (Schumaker SplinePAK)
%tic
[x,y,v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir,...
  dof,A] =  mdsps(xo,yo,v1o,v2o,v3o,e1o,e2o,e3o,...
  ie1o,ie2o,trilo,triro,adjstarto,eadjo,bdyo);
%toc

% refine triangulation (Powell-Sabin-refinement)
TRI = [v1,v2,v3]; 
triplot(TRI,x,y);

%degree of the spline-function
d = 2;

%Regularisationparameter
mu=10;

% 4 options to compute the LSQ-Spline
% -> using Schumaker-Software
% -> using own implementation of the schumaker algorithm (normal equations)
% -> using backslash operator (if mu = 0)
% -> using SVD (prefer this option)

% %%
% tic
% % PLSQ (Schumaker SplinePAK)
% [c_schu,M,t1,t2] = penlsqbiv(d,x,y,v1,v2,v3,...
%     e1,e2,e3,ie1,A,xd,yd,zd,wd,mu);
% t_schu=toc
% 
% % compute energy and RMS-Error of Spline-Functionen
% [E_schu]=energy(x,y,c_schu,v1,v2,v3,e1,e2,e3);
% E_schu=sqrt(E_schu);
% 
% [Nres,maxres_schu]=residuum(TRI,x,y,xd,yd,zd,c_schu,v1,v2,v3,e1,e2,e3);
% rms_schu=Nres/sqrt(length(zd));

%%
tic
%own implementation using normal equations  (PLSQ) für d=2
 [G1,d1,M1,H1]=build_H_d(A,x,y,v1,v2,v3,e1,e2,e3,xd,yd,zd,wd,mu);
 chelp=H1\d1;
 c_ng=A*chelp;
t_ng=toc

% compute energy and RMS-Error of Spline-Functionen
[E_ng]=energy(x,y,c_ng,v1,v2,v3,e1,e2,e3);
E_ng=sqrt(E_ng);

[Nres1,maxres_ng]=residuum(TRI,x,y,xd,yd,zd,c_ng,v1,v2,v3,e1,e2,e3);
rms_ng=Nres1/sqrt(length(zd));

%%
%LSQ-Spline using backslash (just for mu = 0)
if mu==0
%Backslash-Operator
tic
%build observation matrix and compute coefficient vector c
O = build_O(dof,A,x,y,v1,v2,v3,e1,e2,e3,xo,yo,TRIo,xd,yd,ie1o,ie2o);

%compute c (without SVD)
chelp2= O\zd;
c_bs=A*chelp2;
t_bs=toc

% compute energy and RMS-Error of Spline-Functionen
[E_bs]=energy(x,y,c_bs,v1,v2,v3,e1,e2,e3);
E_bs=sqrt(E_bs);

[Nres2,maxres_bs]=residuum(TRI,x,y,xd,yd,zd,c_bs,v1,v2,v3,e1,e2,e3);
rms_bs=Nres2/sqrt(length(zd));
end

%%
%LSQ-Spline using SVD
%tolerance for SVD
tol=10^-6;

tic
%build observation matrix and compute coefficient vector c
Osvd = build_O(dof,A,x,y,v1,v2,v3,e1,e2,e3,xo,yo,TRIo,xd,yd,ie1o,ie2o);

if mu>0
%compute square root of energy matrix -> M^(1/2)
R=sqrtm(M1);
%use SVD on extended Matrix \hat(O) -> here OPLSQ
OPLSQ=[Osvd;sqrt(mu)*R];
[U1,S1,V1]=svd(OPLSQ,'econ');

%determine all singular values (SV) > tol
S1help=1/S1(1,1)*S1;
S1help=diag(S1help);
indsv=find(S1help>tol);

%reduce computation on these SV
S1neu=S1(1:length(indsv),1:length(indsv));
U1neu=U1(:,indsv);
V1neu=V1(:,indsv);

%compute c
chelp3=V1neu/S1neu;
chelp3=chelp3*U1neu'*[zd;zeros(size(R,1),1)];
chelp3=real(chelp3);
csvd=A*chelp3; 

% if mu == 0 -> just using SVD of O
else
[U,S,V]=svd(Osvd,'econ');
chelp3=V/S;
chelp3=chelp3*U'*zd;
csvd=A*chelp3; 
end
tsvd=toc

%energy and rms-error
[Esvd]=energy(x,y,csvd,v1,v2,v3,e1,e2,e3);
Esvd=sqrt(Esvd);

[Nres3,maxres_svd]=residuum(TRI,x,y,xd,yd,zd,csvd,v1,v2,v3,e1,e2,e3);
rms_svd=Nres3/sqrt(length(zd));

%%
% Render the spline (Schumaker SplinePAK) on the rectangle [xmin,xmax]x[ymin,ymax] 
% evaluate on uniform ng x ng grid
ng = 51; 
xmin = min(x)+epsilon*10^k; 
xmax = max(x)-epsilon*10^k; 
ymin = min(y)+epsilon*10^k; 
ymax = max(y)-epsilon*10^k;

% %Schumaker Software
% [xg,yg,g] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c_schu,ng,xmin,xmax,ymin,ymax);
% figure; surfl(xg,yg,g');  colormap(copper); 
% if set==2
%     axis equal;
% end
% title(['Kleinste-Quadrate-Spline mit Strafterm (\mu=',num2str(mu),')']);

% normal equations
[xg1,yg1,g1] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c_ng,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg1,yg1,g1');  colormap(copper); 
if set==2
    axis equal;
end
title(['LSQ-Spline mit Strafterm (eigener Algorithmus, \mu=',num2str(mu),')']);


% SVD
[xg3,yg3,g3] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,csvd,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg3,yg3,g3');  colormap(copper); 
axis equal;
    if set==2
     axis equal;
    end
title('LSQ-Spline (SVD)');

%for mu = 0 -> backslash-Operator
if mu ==0
[xg2,yg2,g2] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c_bs,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg2,yg2,g2');  colormap(copper); 
    if set==2
        axis equal;
    end
title('LSQ-Spline (backslash Operator)');
end

%write information in geo-File
getgeo(csvd,TRI,x,y,v1,v2,v3,e1,e2,e3,ie1);