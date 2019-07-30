% function determines 6 indizes of the mth triangle for getting the
% B-coeffis and the right rows of the Trafomatrix A (just for PS-Space)

%M. Kloppe, Juni 2019

function [index]=finde_ind(m,nv,v1,v2,v3,e1,e2,e3)
% m - Index of the current triangle (of TRIo)
% nv - number of vertices
%v1,v2,v3, e1,e2,e3 - Vectors of indizes of vertices and edges (each
%                     triangle)

% Vector c has for d=2 the following order: 
% c=[coefficients on the vertices; coefficients on the edges]'
%choose right indizes and sort them lexicographically
index=[v1(m);e1(m)+nv;e3(m)+nv;v2(m);e2(m)+nv;v3(m)];
end