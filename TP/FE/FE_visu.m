function FE_visu(UU, mesh, titre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE_visu:
%          
% INPUT * UU (vecteur NbNodes x 1): vecteur solution 
%       * mesh (structure) cf. MESH_build_cartesian
%       * titre (optionel) un titre (string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% le titre est optionnel
if (nargin<3)
    titre = ''; 
end

figure;
trisurf(mesh.AllTriangles,mesh.AllNodes(:,1),mesh.AllNodes(:,2),UU);
view(2);
shading interp
colorbar;

% ajout titre
title(titre, 'Fontsize', 16);
end

