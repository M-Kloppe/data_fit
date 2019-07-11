% Skript zur Ermittlung von Ausgleichsflächen für verstreute Daten

% M. Kloppe, Juni 2019

%Ablauf: (1) Einlesen von Geodaten (Beispiel Stromboli/Schwarzwald)
%        (2) Wahl einer Triangulierung des zugrundeliegenden Gebiets
%        (3) Berechnung von Triangulierungsinformationen sowie MDS und
%            Transformationsmatrix für Powell-Sabin-Raum unter Verwendung
%            von Funktionen aus L. Schumakers SplinePAK (2014)
%            -> SplinePAK: Copyright Larry Schumaker 2014
%        (4) Berechnung des Koeffizientenvektors des gesuchten Splines (LSQ
%            bzw. PLSQ)
%        (5) Plotten der berechneten Spline-Funktion (Schumaker-Software)
%        (6) Berechnung quantitativer Fehlermaße


clear all;
%addpath('splinepak');

%Wähle ein Datenset zum Einlesen der zu verwendeten Daten
%set = 1 -> Stromboli-Daten
%set = 2 -> Schwarzwald-Daten
set=1;

%Einlesen der Geodaten und TRansformation in die richtige Form
%xd,yd,zd - Koordinaten der Geodaten als Spaltenvektoren
%nd -  Anzahl der Datenpunkte
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

%Zugrundeliegendes Gebiet (Rechteck)
epsilon=eps;
k=4;
x1=min(xd)-epsilon*10^k;
y1=min(yd)-epsilon*10^k;
x2=max(xd)+epsilon*10^k;
y2=max(yd)+epsilon*10^k;

%Auswahl der Triangulierung (PDE-Toolbox)
% w=1 grobe Triangulierung (48 Dreiecke)
% w=2 Triangulierung (130 Dreiecke)
% w=3 feine Triangulierung (330 Dreiecke)
w=3;

if w==1
%Triangulierung (mit PDE Toolbox erstellt, 48 Dreiecke) einlesen
Tdat=load('PDEtri1.mat');
TRIo=Tdat.TRIo;
xo = Tdat.xo;
yo = Tdat.yo;
no=length(xo);

%Transformation auf leicht vergrößertes Rechteck
xo=(xo-min(xo))/(max(xo)-min(xo))*abs(x2-x1)+x1;
yo=(yo-min(yo))/(max(yo)-min(yo))*abs(y2-y1)+y1;

elseif w==2  
%Triangulierung (mit PDE Toolbox erstellt, 130 Dreiecke) einlesen
Tdat=load('PDEtri2.mat');
TRIo=Tdat.TRIo;
xo = Tdat.xo;
yo = Tdat.yo;
no=length(xo);

%Transformation auf leicht vergrößertes Rechteck
xo=(xo-min(xo))/(max(xo)-min(xo))*abs(x2-x1)+x1;
yo=(yo-min(yo))/(max(yo)-min(yo))*abs(y2-y1)+y1;

elseif w==3
%Triangulierung (mit PDE Toolbox erstellt, 330 Dreiecke) einlesen
Tdat=load('PDEtri3.mat');
TRIo=Tdat.TRIo;
xo = Tdat.xo;
yo = Tdat.yo;
no=length(xo);

%Transformation auf leicht vergrößertes Rechteck
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

% Simulation von Rauscheffekten
e = 0;
if e > 0
 %erzeuge Gaußsches Rauschen
 %rv=-1+2.*rand(nd,1);
 %zd = zd + e*rv;
 
 %Lade bereits erzeugten Rauschvektor
 noise=load('rauschen2.mat');
 zd = zd + e*noise.rv;
end

%Find a minimal determing set and transformation matrix (Schumaker SplinePAK)
%tic
[x,y,v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir,...
  dof,A] =  mdsps(xo,yo,v1o,v2o,v3o,e1o,e2o,e3o,...
  ie1o,ie2o,trilo,triro,adjstarto,eadjo,bdyo);
%toc

% Verfeinerte Triangulierung (Powell-Sabin-Verfeinerung)
TRI = [v1,v2,v3]; 

%Grad der Spline-Funktionen
d = 2;

%Regularisierungsparameter
mu=1000;

%%
tic
% PLSQ (Schumaker SplinePAK)
[c_schu,M,t1,t2] = penlsqbiv(d,x,y,v1,v2,v3,...
    e1,e2,e3,ie1,A,xd,yd,zd,wd,mu);
t_schu=toc

% berechne Energie und RMS-Error der Spline-Funktionen (->schumaker)
[E_schu]=energy(x,y,c_schu,v1,v2,v3,e1,e2,e3);
E_schu=sqrt(E_schu);

[Nres,maxres_schu]=residuum(TRI,x,y,xd,yd,zd,c_schu,v1,v2,v3,e1,e2,e3);
rms_schu=Nres/sqrt(length(zd));

%%
tic
%eigener Algorithmus (PLSQ) für d=2
 [G1,d1,M1,H1]=build_H_d(A,x,y,v1,v2,v3,e1,e2,e3,xd,yd,zd,wd,mu);
 chelp=H1\d1;
 c_ng=A*chelp;
t_ng=toc

% berechne Energie und RMS-Error der Spline-Funktionen (->Normalgl.)
[E_ng]=energy(x,y,c_ng,v1,v2,v3,e1,e2,e3);
E_ng=sqrt(E_ng);

[Nres1,maxres_ng]=residuum(TRI,x,y,xd,yd,zd,c_ng,v1,v2,v3,e1,e2,e3);
rms_ng=Nres1/sqrt(length(zd));

%%
%LSQ-Spline mit backslash (nur bei mu = 0 ausführen)
if mu==0
%Backslash-Operator
tic
%Aufstellen Beobachtungsmatrix und Bestimmung c
O = build_O(dof,A,x,y,v1,v2,v3,e1,e2,e3,xo,yo,TRIo,xd,yd,ie1o,ie2o);

%Berechnung von c (Ohne SVD)
chelp2= O\zd;
c_bs=A*chelp2;
t_bs=toc

%Energie und rms
[E_bs]=energy(x,y,c_bs,v1,v2,v3,e1,e2,e3);
E_bs=sqrt(E_bs);

[Nres2,maxres_bs]=residuum(TRI,x,y,xd,yd,zd,c_bs,v1,v2,v3,e1,e2,e3);
rms_bs=Nres2/sqrt(length(zd));
end

%%
%LSQ-Spline mit SVD
%Toleranznieveau fuer SVD
tol=10^-6;

tic
%Aufstellen Beobachtungsmatrix und Bestimmung c
Osvd = build_O(dof,A,x,y,v1,v2,v3,e1,e2,e3,xo,yo,TRIo,xd,yd,ie1o,ie2o);

if mu>0
%Berechnung M^(1/2)
R=sqrtm(M1);
%Löse über SVD der Erweiterten Matrix \hat(O) -> hier OPLSQ
OPLSQ=[Osvd;sqrt(mu)*R];
[U1,S1,V1]=svd(OPLSQ,'econ');

%ermittle alle Singulärwerte (SW) > tol
S1help=1/S1(1,1)*S1;
S1help=diag(S1help);
indsv=find(S1help>tol);

%Reduziere Berechnungsformel auf diese SW
S1neu=S1(1:length(indsv),1:length(indsv));
U1neu=U1(:,indsv);
V1neu=V1(:,indsv);

%Berechne c
chelp3=V1neu/S1neu;
chelp3=chelp3*U1neu'*[zd;zeros(size(R,1),1)];
chelp3=real(chelp3);
csvd=A*chelp3; 

else
%Löse über SVD der Beobachtungsmatrix O
[U,S,V]=svd(Osvd,'econ');
chelp3=V/S;
chelp3=chelp3*U'*zd;
csvd=A*chelp3; 
end
tsvd=toc

%Energie und rms bei SVD
[Esvd]=energy(x,y,csvd,v1,v2,v3,e1,e2,e3);
Esvd=sqrt(Esvd);

[Nres3,maxres_svd]=residuum(TRI,x,y,xd,yd,zd,csvd,v1,v2,v3,e1,e2,e3);
rms_svd=Nres3/sqrt(length(zd));

%%
% Render the spline (Schumaker SplinePAK)
% müssen auch hier Rechteck anpassen [xmin,xmax]x[ymin,ymax] 
% werte Spline über uniformen Gitter mit ng x ng Punkten aus
ng = 51; 
xmin = min(x)+epsilon*10^k; 
xmax = max(x)-epsilon*10^k; 
ymin = min(y)+epsilon*10^k; 
ymax = max(y)-epsilon*10^k;

% mit Schumaker Software berechnet
[xg,yg,g] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c_schu,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg,yg,g');  colormap(copper); 
if set==2
    axis equal;
end
title(['Kleinste-Quadrate-Spline mit Strafterm (\mu=',num2str(mu),')']);

% mit eigener Routine berechnet - Normalgleichungen
[xg1,yg1,g1] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c_ng,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg1,yg1,g1');  colormap(copper); 
if set==2
    axis equal;
end
title(['LSQ-Spline mit Strafterm (eigener Algorithmus, \mu=',num2str(mu),')']);


% mit eigener Routine berechnet - SVD
[xg3,yg3,g3] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,csvd,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg3,yg3,g3');  colormap(copper); 
if set==2
    axis equal;
end
title('LSQ-Spline (SVD)');

%für mu = 0 plotte auch Spline bei Verwendung des backslash-Operators
if mu ==0
% mit eigener Routine berechnet - Backslash
[xg2,yg2,g2] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c_bs,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg2,yg2,g2');  colormap(copper); 
if set==2
    axis equal;
end
title('LSQ-Spline (backslash Operator)');
end
