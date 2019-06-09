% Skript zur Ermittlung von Ausgleichsfl�chen f�r verstreute Daten

% M. Kloppe, Juni 2019

%Ablauf: (1) Einlesen von Geodaten (Beispiel Stromboli/Schwarzwald)
%        (2) Wahl einer Triangulierung des zugrundeliegenden Gebiets
%        (3) Berechnung von Triangulierungsinformationen sowie MDS und
%            Transformationsmatrix f�r Powell-Sabin-Raum unter Verwendung
%            von Funktionen aus L. Schumakers SplinePAK (2014)
%            -> SplinePAK: Copyright Larry Schumaker 2014
%        (4) Berechnung des Koeffizientenvektors des gesuchten Splines (LSQ
%            bzw. PLSQ)
%        (5) Plotten der berechneten Spline-Funktion (Schumaker-Software)
%        (6) Berechnung quantitativer Fehlerma�e


clear all;
addpath('splinepak');

%W�hle ein Datenset zum Einlesen der zu verwendeten Daten
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

%Transformation auf leicht vergr��ertes Rechteck
xo=(xo-min(xo))/(max(xo)-min(xo))*abs(x2-x1)+x1;
yo=(yo-min(yo))/(max(yo)-min(yo))*abs(y2-y1)+y1;

elseif w==2  
%Triangulierung (mit PDE Toolbox erstellt, 130 Dreiecke) einlesen
Tdat=load('PDEtri2.mat');
TRIo=Tdat.TRIo;
xo = Tdat.xo;
yo = Tdat.yo;
no=length(xo);

%Transformation auf leicht vergr��ertes Rechteck
xo=(xo-min(xo))/(max(xo)-min(xo))*abs(x2-x1)+x1;
yo=(yo-min(yo))/(max(yo)-min(yo))*abs(y2-y1)+y1;

elseif w==3
%Triangulierung (mit PDE Toolbox erstellt, 330 Dreiecke) einlesen
Tdat=load('PDEtri3.mat');
TRIo=Tdat.TRIo;
xo = Tdat.xo;
yo = Tdat.yo;
no=length(xo);

%Transformation auf leicht vergr��ertes Rechteck
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
 %erzeuge Gau�sches Rauschen
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
mu=0;

tic
% PLSQ (Schumaker SplinePAK)
[c,M,t1,t2] = penlsqbiv(d,x,y,v1,v2,v3,...
    e1,e2,e3,ie1,A,xd,yd,zd,wd,mu);
t=toc

tic
%eigener Algorithmus (PLSQ) f�r d=2
 [G1,d1,M1,H1]=build_H_d(A,x,y,v1,v2,v3,e1,e2,e3,xd,yd,zd,wd,mu);
 chelp=H1\d1;
 c1=A*chelp;
ts=toc

%Alternativen f�r mu = 0 -> backslash und SVD (nur bei mu = 0 ausf�hren)
if mu==0
%Backslash-Operator
tic
%Aufstellen Beobachtungsmatrix und Bestimmung c
O = build_O(dof,A,x,y,v1,v2,v3,e1,e2,e3,xo,yo,TRIo,xd,yd,ie1o,ie2o);

%Berechnung von c (Ohne SVD)
chelp2= O\zd;
c2=A*chelp2;
tO=toc

%Energie und rms
[E2]=energy(x,y,c2,v1,v2,v3,e1,e2,e3);
E2=sqrt(E2);

[Nres2,maxres2]=residuum(TRI,x,y,xd,yd,zd,c2,v1,v2,v3,e1,e2,e3);
rms2=Nres2/sqrt(length(zd));

%SVD
tic
%Aufstellen Beobachtungsmatrix und Bestimmung c
Osvd = build_O(dof,A,x,y,v1,v2,v3,e1,e2,e3,xo,yo,TRIo,xd,yd,ie1o,ie2o);

[U,S,V]=svd(Osvd,'econ');
r=size(S,1);
chelp3=zeros(size(V,1),1);
for i=1:r
    chelp3=chelp3+1/S(i,i)*U(:,i)'*zd*V(:,i);
end
c3=A*chelp3; 
tsvd=toc

%Energie und rms
[E3]=energy(x,y,c3,v1,v2,v3,e1,e2,e3);
E3=sqrt(E3);

[Nres3,maxres3]=residuum(TRI,x,y,xd,yd,zd,c3,v1,v2,v3,e1,e2,e3);
rms3=Nres3/sqrt(length(zd));
end


% Render the spline (Schumaker SplinePAK)
% m�ssen auch hier Rechteck anpassen [xmin,xmax]x[ymin,ymax] 
% werte Spline �ber uniformen Gitter mit ng x ng Punkten aus
ng = 51; 
xmin = min(x)+epsilon*10^k; 
xmax = max(x)-epsilon*10^k; 
ymin = min(y)+epsilon*10^k; 
ymax = max(y)-epsilon*10^k;

% mit Schumaker Software berechnet
[xg,yg,g] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg,yg,g');  colormap(copper); 
if set==2
    axis equal;
end
title(['Kleinste-Quadrate-Spline mit Strafterm (\mu=',num2str(mu),')']);

% mit eigener Routine berechnet
[xg1,yg1,g1] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c1,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg1,yg1,g1');  colormap(copper); 
if set==2
    axis equal;
end
title(['LSQ-Spline mit Strafterm (eigener Algorithmus, \mu=',num2str(mu),')']);

% berechne Energie der Spline-Funktionen
[E]=energy(x,y,c,v1,v2,v3,e1,e2,e3);
[E1]=energy(x,y,c1,v1,v2,v3,e1,e2,e3);

E=sqrt(E);
E1=sqrt(E1);

%Berechne RMS-Error (=Euklid. Norm des Residuums/sqrt(Anzahl Datenpunkte)) 
%und maximale Abweichung 
[Nres,maxres]=residuum(TRI,x,y,xd,yd,zd,c,v1,v2,v3,e1,e2,e3);
rms=Nres/sqrt(length(zd));

[Nres1,maxres1]=residuum(TRI,x,y,xd,yd,zd,c1,v1,v2,v3,e1,e2,e3);
rms1=Nres1/sqrt(length(zd));

%f�r mu = 0 plotte auch Splines der anderen zwei Varianten
if mu ==0
% mit eigener Routine berechnet - Backslash
[xg2,yg2,g2] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c2,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg2,yg2,g2');  colormap(copper); 
if set==2
    axis equal;
end
title('LSQ-Spline (backslash Operator)');

% mit eigener Routine berechnet - SVD
[xg3,yg3,g3] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c3,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg3,yg3,g3');  colormap(copper); 
if set==2
    axis equal;
end
title('LSQ-Spline (SVD)');
end

