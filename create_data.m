%function changes the shape of a given matrix of data (Matrix B - Stromboli)
%-> determines vectors xd,yd,zd with the cartesian coords of the data points 

%M. Kloppe, Juni 2019
function [xd,yd,zd,nd]=create_data(B)
%Dimension of matrix B
[z,s]=size(B);

%in matrix B -> points are something like pixels in xy-plane
% create xd and yd using linspace
xhelp=linspace(1,s,s);
yhelp=linspace(1,z,z);

%build uniform grid and reshape elements to get vectors
[X,Y]=meshgrid(xhelp,yhelp);
xd=reshape(X,[length(xhelp)*length(yhelp),1]);
yd=reshape(Y,[length(xhelp)*length(yhelp),1]);

%build vector of the geographical heights
zd=reshape(B,[length(xhelp)*length(yhelp),1]);
zd=double(zd);
nd=length(zd);
end