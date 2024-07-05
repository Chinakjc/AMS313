function Mel = FE_matM_elem(node1_coord, node2_coord, node3_coord)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul la matrice de masse elementaire en P1 lagrange
%
% INPUT * node1_coord, node1_coord, node1_coord (chacun de taille 2) 
%         les coordonnees (x,y) des 3 noeuds formant le triangle
%
% OUTPUT - Mel (taille 3x3) matrice de masse elementaire 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordonnees des 3 noeuds formant le triangle
x1 = node1_coord(1); y1 = node1_coord(2);
x2 = node2_coord(1); y2 = node2_coord(2);
x3 = node3_coord(1); y3 = node3_coord(2);

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

% calcul de la matrice de masse
% -----------------------------
Mel = zeros(3,3);
for i=1:3
    for j=1:3
        if i==j
            Mel(i,j) = 2/24*abs(D);
        else
            Mel(i,j) = 1/24*abs(D);
        end
    end
end

end
