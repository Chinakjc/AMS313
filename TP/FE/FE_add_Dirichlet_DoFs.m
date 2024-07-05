function [ UU_full ] = FE_add_Dirichlet_DoFs( UU, mesh , DofNodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE_relevement:
%          
% INPUT * UU (taille NbDof) : solution sans les noeuds 
%                             de Dirichlet homogene
%       * mesh (structure) : le maillage
%       * DofNodes (taille NbDof) : indices des noeuds correspondant a des 
%                 degres de liberte du probleme
%
% OUTPUT * UU_full (taille NbNodes) : solution avec les noeuds de Dirichlet
%                                     homogene
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NbDof = size(DofNodes,1);
if (size(UU,1) ~= NbDof)
    error('le vecteur en entree pas de la bonne taille');
end

NbNodes = mesh.NbNodes;
UU_full = zeros(NbNodes,1);
% UU_full et UU coincident sur les degres de liberte du probleme,
% UU_full vaut 0 en dehors des degres de liberte du probleme
UU_full(DofNodes) = UU;

end

