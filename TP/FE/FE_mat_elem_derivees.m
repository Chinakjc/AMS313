function [DXel, DYel] =FE_mat_elem_derivees(node1_coord, node2_coord, node3_coord)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul les matrices elementaires pour les derivees
%          
% INPUT * node1_coord, node1_coord, node1_coord (chacun de taille 2) 
%         les coordonnees (x,y) des 3 noeuds formant le triangle
%
% OUTPUT - DXel, DYel matrices elementaires 3x3 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordonnees des 3 noeuds
x1 = node1_coord(1); y1 = node1_coord(2);
x2 = node2_coord(1); y2 = node2_coord(2);
x3 = node3_coord(1); y3 = node3_coord(2);

% Determinant
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
% signe du determinant
sgnD = abs(D) / D;

% Declaration et remplisage des matrice elementaires.
%---------------------------------------------------
DXel = zeros(3,3);
DYel = zeros(3,3);

% premiere colonne DXel
DXel(1,1)=(y2-y3) * sgnD/6.0;
DXel(2,1)=(y2-y3) * sgnD/6.0;
DXel(3,1)=(y2-y3) * sgnD/6.0;
% deuxieme colonne DXel
DXel(1,2)=(y3-y1) * sgnD/6.0;
DXel(2,2)=(y3-y1) * sgnD/6.0;
DXel(3,2)=(y3-y1) * sgnD/6.0;
% troisieme colonne DXel
DXel(1,3)=(y1-y2) * sgnD/6.0;
DXel(2,3)=(y1-y2) * sgnD/6.0;
DXel(3,3)=(y1-y2) * sgnD/6.0;

% premiere colonne DYel
DYel(1,1)=(x3-x2) * sgnD/6.0;
DYel(2,1)=(x3-x2) * sgnD/6.0;
DYel(3,1)=(x3-x2) * sgnD/6.0;
% deuxieme colonne DYel
DYel(1,2)=(x1-x3) * sgnD/6.0;
DYel(2,2)=(x1-x3) * sgnD/6.0;
DYel(3,2)=(x1-x3) * sgnD/6.0;
% troisieme colonne DYel
DYel(1,3)=(x2-x1) * sgnD/6.0;
DYel(2,3)=(x2-x1) * sgnD/6.0;
DYel(3,3)=(x2-x1) * sgnD/6.0;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
