function [Kel] = FE_matK_elem(node1_coord, node2_coord, node3_coord)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul la matrice de raideur elementaire en P1 lagrange
%          
% INPUT * node1_coord, node1_coord, node1_coord (chacun de taille 2) 
%         les coordonnees (x,y) des 3 noeuds formant le triangle
%
% OUTPUT - Kel (taille 3x3) matrice de raideur elementaire 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordonnees des 3 noeuds formant le triangle
x1 = node1_coord(1); y1 = node1_coord(2);
x2 = node2_coord(1); y2 = node2_coord(2);
x3 = node3_coord(1); y3 = node3_coord(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;


% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
for i=1:3
  for j=1:3
    Kel(i,j) = dot(norm(i,:),norm(j,:))/(2*abs(D));
  end
end

end
