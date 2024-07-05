function [S1el, S2el, S3el, S4el, S5el, S6el] = FE_mat_elem_surface(node1_coord, node2_coord, ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul les matrices de masse surfaciques pour les frontieres 
%
% INPUT * node1_coord : coordonnees x,y du noeud 
%       * node2_coord : coordonnees x,y du noeud 
%       * ref: reference de l'arete
%
% OUTPUT - S1el matrice de masse surfacique elementaire pour Gamma1 (matrice 2x2)
%		 - S2el matrice de masse surfacique elementaire pour Gamma2 (matrice 2x2)
%        - S3el matrice de masse surfacique elementaire pour Gamma3 (matrice 2x2)
%        - S4el matrice de masse surfacique elementaire pour Gamma4 (matrice 2x2)
%        - S5el matrice de masse surfacique elementaire pour Gamma5 (matrice 2x2)
%        - S6el matrice de masse surfacique elementaire pour Gamma6 (matrice 2x2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordonnees des 2 noeuds formant l'arete
x1 = node1_coord(1); y1 = node1_coord(2);
x2 = node2_coord(1); y2 = node2_coord(2);

% Longueur de l'arete.
Long = sqrt((x2-x1)^2+(y2-y1)^2);

% Declarations
S1el = zeros(2,2);
S2el = zeros(2,2);
S3el = zeros(2,2);
S4el = zeros(2,2);
S5el = zeros(2,2);
S6el = zeros(2,2);
% caclul
if ref==1      % Gamma1 
    S1el = [Long/3 Long/6; Long/6 Long/3];
elseif ref==2  % Gamma2 
    S2el = [Long/3 Long/6; Long/6 Long/3];
elseif ref==3  % Gamma3
    S3el = [Long/3 Long/6; Long/6 Long/3];  
elseif ref==4  % Gamma4 
    S4el = [Long/3 Long/6; Long/6 Long/3];  
elseif ref==5  % Gamma5 
    S5el = [Long/3 Long/6; Long/6 Long/3];
elseif ref==6  % Gamma6 
    S6el = [Long/3 Long/6; Long/6 Long/3]; 
else
    error('la reference nest pas definie')
end

end