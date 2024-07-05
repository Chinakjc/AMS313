function MESH_visu(mesh, titre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualisation d'un maillage triangulaire 2D
%
% INPUT  * mesh (structure) : objet maillage
%        * titre (optionel) un titre (string)
%
% OUTPUT une fenetre graphique
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control on the input args
if (nargin<2), titre = ''; end;

%visualisation du maillage
figure;
hold on

% maillage
NbNodes = mesh.NbNodes;
trimesh(mesh.AllTriangles,mesh.AllNodes(:,1),mesh.AllNodes(:,2),zeros(NbNodes,1));
view(2);
axis('equal');

% ajout du titre
title(titre);

hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
