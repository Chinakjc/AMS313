function FE_visu_vector( VX,VY, mesh, titre )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE_visu_vector:
%          
% INPUT * VX (vecteur NbNodes x 1): vecteur des composantes X 
%       * VY (vecteur NbNodes x 1): vecteur des composantes Y 
%       * mesh (structure) cf. MESH_build_cartesian
%       * titre (optionel) un titre (string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% le titre est optionnel
if (nargin<4)
    titre = ''; 
end

figure;
NbNodes = mesh.NbNodes;
nb_arrows = 0;
restriction = [];
if (NbNodes < 500)
    nb_arrows = NbNodes;
    restriction = 1:nb_arrows;
else
    nb_arrows = 500;
    restriction = floor(1 + (NbNodes - 1)*rand(nb_arrows,1));
end
quiver(mesh.AllNodes(restriction,1),  mesh.AllNodes(restriction,2), ...
       VX(restriction), VY(restriction));

% ajout titre
title(titre, 'Fontsize', 16);
end

