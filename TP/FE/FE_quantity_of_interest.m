function [ JJ ] = FE_quantity_of_interest( mesh, DofNodes )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE_assemblages :
% Assemblage des matrices elements finis.
% Assemblage du second membre elements finis.
%          
% INPUT * mesh (structure) : tous les tableaux concernant le maillage
%                            (cf. MESH_build_cartesian.m)
%       * DofNodes (taille NbDof) : indices des noeuds qui correspondent a
%              des degres de liberte (= tous les noeuds a l'exception des 
%              noeuds de dirichlet)
%
% OUTPUT - JJ (taille NbDof) : forme lineaire qui integre la fonction sur 
%                  le bord Gamma6
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extraction des tableaux utiles depuis la structure mesh
NbNodes= mesh.NbNodes;
AllNodes = mesh.AllNodes;
NbBoundaryEdges = mesh.NbBoundaryEdges;
AllBoundaryEdges= mesh.AllBoundaryEdges;
RefBoundaryEdges= mesh.RefBoundaryEdges;

% matrice de masse surfacique sur Gamma1 U Gamma2.
SS_in = sparse(NbNodes,NbNodes); 
% matrice de masse surfacique sur Gamma4 U Gamma5.
SS_out = sparse(NbNodes,NbNodes); 
% assemblage
for a=1:NbBoundaryEdges
    node1_coords = AllNodes(AllBoundaryEdges(a,1),:);
    node2_coords = AllNodes(AllBoundaryEdges(a,2),:);
    % Calcul des matrices de masse sur la surface Gamma1 U Gamma2
    Sinel = zeros(2,2);
    ref_a = RefBoundaryEdges(a,1);
    if (ref_a==1 || ref_a==2 )
        % coordonnees des 2 noeuds formant l'arete
        x1 = node1_coords(1); y1 = node1_coords(2);
        x2 = node2_coords(1); y2 = node2_coords(2);
        % Longueur de l'arete.
        Long = sqrt((x2-x1)^2+(y2-y1)^2);
        % matrice elem
        Sinel = [Long/3 Long/6; Long/6 Long/3]; 
    end
    % Calcul des matrices de masse sur la surface Gamma5 U Gamma4
    Soutel = zeros(2,2);
    if (ref_a==4 || ref_a==5 )
        % coordonnees des 2 noeuds formant l'arete
        x1 = node1_coords(1); y1 = node1_coords(2);
        x2 = node2_coords(1); y2 = node2_coords(2);
        % Longueur de l'arete.
        Long = sqrt((x2-x1)^2+(y2-y1)^2);
        % matrice elem
        Soutel = [Long/3 Long/6; Long/6 Long/3]; 
    end
    % On fait l'assemblage de la matrice globale SS.
    for i=1:2
        I = AllBoundaryEdges(a,i);
        for j=1:2
            J = AllBoundaryEdges(a,j);
            SS_in(I,J) = SS_in(I,J)   + Sinel(i,j);
            SS_out(I,J) = SS_out(I,J) + Soutel(i,j);
        end 
    end 
end 

% formes lineaires
% --------------
% forme qui integre sur la surface Gamma1 U Gamma2
JJ_in = SS_in*ones(NbNodes,1);
% forme qui integre sur la surface Gamma4 U Gamma5
JJ_out = SS_out*ones(NbNodes,1);

% quantite d'interet = (integrale sur Gamma4 U Gamma5) - (integrale sur Gamma1 U Gamma2) 
JJ = JJ_out - JJ_in;
% extraction des noeuds de Dirichlet
JJ = JJ(DofNodes);

end

