function PHIr = FE_relevement( mesh )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE_relevement:
%          
% INPUT * mesh (structure) : le maillage
%
% OUTPUT * PHIr (taille NbNodes) : relevement pour la cond. de Dirichlet
%                   non-homogene gdir
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PHIr = zeros(mesh.NbNodes,1);
% construction du relevement
for i=1:mesh.NbNodes
    xi=mesh.AllNodes(i,1);
    yi=mesh.AllNodes(i,2);
    PHIr(i,1) = FE_gdir(xi,yi);
end

end

