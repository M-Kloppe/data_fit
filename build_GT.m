%function computes 6x6 Matrix G(T) and vector d(T) for algorithm 7.3
%(Schumaker Computational Methods) for d=2

%M. Kloppe, Juni 2019

%TRI ... triangulation (triangulation object!)
%index ... index of triangle T in TRI
%wd ... weights
%xd,yd,zd ... data points (in columns)

function [GT,dT]=build_GT(TRI,index,wd,xd,yd,zd)
%Initialize
GT=zeros(6,6);
dT=zeros(6,1);

%Index as a vector of length x -> use for cartesianToBarycentric
ihelp=index*ones(length(xd),1);

%transform data points in barycentric coordinates
B = cartesianToBarycentric(TRI,ihelp,[xd,yd]);

%use quadratic interpolation
d=2;

%indices and normalisation factor of Bernstein-Polynomials
d_f=factorial(d);
M=[2 0 0;0 2 0;0 0 2;1 1 0;1 0 1;0 1 1];
D=factorial(M);

%Bernstein-Basispolynomials (lexikographical order)
bern1=@(b1,b2,b3) d_f/(D(1,1)*D(1,2)*D(1,3)).*b1.^M(1,1).*b2.^M(1,2).*b3.^M(1,3);
bern4=@(b1,b2,b3) d_f/(D(2,1)*D(2,2)*D(2,3)).*b1.^M(2,1).*b2.^M(2,2).*b3.^M(2,3);
bern6=@(b1,b2,b3) d_f/(D(3,1)*D(3,2)*D(3,3)).*b1.^M(3,1).*b2.^M(3,2).*b3.^M(3,3);
bern2=@(b1,b2,b3) d_f/(D(4,1)*D(4,2)*D(4,3)).*b1.^M(4,1).*b2.^M(4,2).*b3.^M(4,3);
bern3=@(b1,b2,b3) d_f/(D(5,1)*D(5,2)*D(5,3)).*b1.^M(5,1).*b2.^M(5,2).*b3.^M(5,3);
bern5=@(b1,b2,b3) d_f/(D(6,1)*D(6,2)*D(6,3)).*b1.^M(6,1).*b2.^M(6,2).*b3.^M(6,3);


%Evaluate Basispolynomials for the data
B1=bern1(B(:,1),B(:,2),B(:,3));
B2=bern2(B(:,1),B(:,2),B(:,3));
B3=bern3(B(:,1),B(:,2),B(:,3));
B4=bern4(B(:,1),B(:,2),B(:,3));
B5=bern5(B(:,1),B(:,2),B(:,3));
B6=bern6(B(:,1),B(:,2),B(:,3));

%all the vectors in on matrix
Bmat=[B1,B2,B3,B4,B5,B6];
W=diag(wd);

%Update GT
GT=Bmat'*W*Bmat;

%compute rT
r1=zd'*diag(wd)*B1;
r2=zd'*diag(wd)*B2;
r3=zd'*diag(wd)*B3;
r4=zd'*diag(wd)*B4;
r5=zd'*diag(wd)*B5;
r6=zd'*diag(wd)*B6;

dT=[r1;r2;r3;r4;r5;r6];
end