% Funktion ermittelt die 6 Indizes des m-ten Dreiecks für die Zuordnung der
% B-Koeffis bzw. der Zeilen von der Trafomatrix A (geeignet bei Verwendung
% des Powell-Sabin-Raumes)

%M. Kloppe, Juni 2019

function [index]=finde_ind(m,nv,v1,v2,v3,e1,e2,e3)
% m - Index des aktuellen Weltdreiecks
% nv - Anzahl aller Eckpunkte
% v1,v2,v3 - Liste der Eckenindizes der Dreiecke
% e1,e2,e3 - Liste der Kantenindizes der Dreiecke

% Vektor c hat für d=2 die Ordnung: 
% c=[alle Koeffis zu den Eckpunkten; alle Koeffis zu Punkten auf Kanten]'
%Wähle richtige Indizes aus und ordne diese lexikographisch an
index=[v1(m);e1(m)+nv;e3(m)+nv;v2(m);e2(m)+nv;v3(m)];
end