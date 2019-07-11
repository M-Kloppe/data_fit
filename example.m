% Failing example Stromboli (very coarse mesh)
f = matfile('geotestdata.mat');
getgeo(f.c, f.TRI, f.x, f.y, f.v1, f.v2, f.v3, f.e1, f.e2, f.e3, f.ie1);

% Better example ?
CP_plsq_ps;
getgeo(csvd, TRI, x, y, v1, v2, v3, e1, e2, e3, ie1);
