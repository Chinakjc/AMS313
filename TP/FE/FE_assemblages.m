function [ DofNodes, AA, LL, ...
           MM, DDX, DDY, BB, AA_decomp, LL_decomp]  = FE_assemblages(mesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE_assemblages :
% Assemblage des matrices elements finis.
% Assemblage du second membre elements finis.
%          
% INPUT * mesh (structure) : tous les tableaux concernant le maillage
%                            (cf. MESH_build_cartesian.m)
%
% OUTPUT - DofNodes (taille NbDof) : indices des noeuds qui correspondent a
%              des degres de liberte (= tous les noeuds a l'exception des 
%              noeuds de dirichlet)
%        - AA (taille NbDof x NbDof) : matrice de l'operateur -div(K Grad) 
%                                      + cond. Dirichlet 
%        - LL (taille NbDof) : second membre 
%        - MM  (taille NbNodes x NbNodes): matrice masse 
%        - DDX (taille NbNodes x NbNodes): matrice derivation selon x
%        - DDY (taille NbNodes x NbNodes): matrice derivation selon y 
%        - BB (taille NbDof x NbDof) : matrice du produit scalaire H1 
%                                           + cond. Dirichlet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('----------------------------');
disp(' Assemblage des matrices EF');
disp('----------------------------');

% recuperation de tous les tableaux a partir de la structure mesh
NbNodes = mesh.NbNodes;            % nombre de noeuds
NbTriangles = mesh.NbTriangles;    % nombre de triangle
AllNodes = mesh.AllNodes;          % (taille NbNodes x 2): coordonnees x,y 
%                                    des noeuds
AllTriangles = mesh.AllTriangles;  % (taille NbTriangles x 3): indices des
%                                    3 noeuds de chaque triangle
RefNodes = mesh.RefNodes;          % (taille NbNodes) : references des noeuds
RefTriangles = mesh.RefTriangles;  % (taille NbTriangles) : references des 
%                                    triangles
NbBoundaryEdges = mesh.NbBoundaryEdges;   % nombre d'aretes de bord
AllBoundaryEdges = mesh.AllBoundaryEdges; % (taille NbBoundaryEdges x 2): 
%                                  indices 2 noeuds  de chaque arete de bord
RefBoundaryEdges = mesh.RefBoundaryEdges; % (taille NbBoundaryEdges): references des aretes
%                                   de bord


disp(strcat('Nombre de noeuds = ', num2str(NbNodes)));

% Construction des matrices volumiques sur les sous-domaines Omega_j
% -----------------------------------------------------------------
% Note: Le maillage est tel que la reference des triangles de Omega_j est j
AA = sparse(NbNodes,NbNodes); % matrice de l'operateur -div(K grad)
KK = sparse(NbNodes,NbNodes); % matrice de rigidite
MM = sparse(NbNodes,NbNodes); % matrice de masse
DDX = sparse(NbNodes,NbNodes); % matrice pour la derivee par rapport a X
DDY = sparse(NbNodes,NbNodes); % matrice pour la derivee par rapport a Y
AA_decomp = cell(3); % decomposition de matrice AA par A0, A1 et A2
LL_decomp = cell(3); % decomposition de seconmembre LL par L0, L1 et L2
AA0 = sparse(NbNodes,NbNodes);% A0
AA1 = sparse(NbNodes,NbNodes);% A1
AA2 = sparse(NbNodes,NbNodes);% A2

% assemblages
for l=1:NbTriangles
  % Coordonnees des sommets du triangle
  node1_coords = AllNodes(AllTriangles(l,1),:);
  node2_coords = AllNodes(AllTriangles(l,2),:);
  node3_coords = AllNodes(AllTriangles(l,3),:);
  % Rigidite elementaire
  Kel = FE_matK_elem(node1_coords, node2_coords, node3_coords);
  % La matrice elementaire de A est egale a la matrice de rigidite
  % elementaire multipliee par le coeff de diffusion kappa
  kappa = FE_kappa_coeff(RefTriangles(l));
  Ael = kappa * Kel;
  % Masse elementaire
  Mel = FE_matM_elem(node1_coords, node2_coords, node3_coords);
  % Derivees /x et /y elementaires
  [DXel, DYel] =FE_mat_elem_derivees(node1_coords, node2_coords, node3_coords);
  % Assemblage des matrices globales 
  for i=1:3
      I = AllTriangles(l,i);
      for j=1:3
          J = AllTriangles(l,j);
          MM(I,J) = MM(I,J) + Mel(i,j);
          KK(I,J) = KK(I,J) + Kel(i,j);
          AA(I,J) = AA(I,J) + Ael(i,j);
          DDX(I,J)= DDX(I,J)+ DXel(i,j);  
          DDY(I,J)= DDY(I,J)+ DYel(i,j);
          switch (RefTriangles(l))
            case 0
                AA0(I,J) = AA0(I,J) + Kel(i,j); % omega 0
            case 1
                AA1(I,J) = AA1(I,J) + Kel(i,j); % omega 1
            case 2
                AA2(I,J) = AA2(I,J) + Kel(i,j); % omega 2
            otherwise
                error('la reference doit etre 0,1 ou 2 !!!')
        end
      end
  end
end 

% matrice du produit scalaire H1
BB = KK + MM;


% Construction des matrices surfaciques sur les frontieres Gamma_j
% -----------------------------------------------------------------
% Note: Le maillage est tel que la reference des aretes de Gamma_j est j
SS1 = sparse(NbNodes,NbNodes); % matrice de masse surfacique sur Gamma1.
SS2 = sparse(NbNodes,NbNodes); % matrice de masse surfacique sur Gamma2.
SS3 = sparse(NbNodes,NbNodes); % matrice de masse surfacique sur Gamma3.
SS4 = sparse(NbNodes,NbNodes); % matrice de masse surfacique sur Gamma4.
SS5 = sparse(NbNodes,NbNodes); % matrice de masse surfacique sur Gamma5.
SS6 = sparse(NbNodes,NbNodes); % matrice de masse surfacique sur Gamma6.
% assemblage
for a=1:NbBoundaryEdges
    node1_coords = AllNodes(AllBoundaryEdges(a,1),:);
    node2_coords = AllNodes(AllBoundaryEdges(a,2),:);
    % Calcul des matrices de masse sur les surfaces Gammaj
    [S1el,S2el,S3el,S4el,S5el,S6el] = ...
        FE_mat_elem_surface(node1_coords, node2_coords,RefBoundaryEdges(a,1));
    % On fait l'assemblage des matrices globale SS1, SS2, ..., SS6.
    for i=1:2
        I = AllBoundaryEdges(a,i);
        for j=1:2
            J = AllBoundaryEdges(a,j);
            SS1(I,J) = SS1(I,J) + S1el(i,j);
            SS2(I,J) = SS2(I,J) + S2el(i,j);
            SS3(I,J) = SS3(I,J) + S3el(i,j);
            SS4(I,J) = SS4(I,J) + S4el(i,j);
            SS5(I,J) = SS5(I,J) + S5el(i,j);
            SS6(I,J) = SS6(I,J) + S6el(i,j);
        end 
    end 
end 

% Assemblage second membre
% ------------------------
% Note: on rappelle le schema du domaine :
%      G2    G3    G4
%    |----|------|----|
% G1 |     \Omega     | G5
%    |----------------|
%            G6
%  - Condition de Neumann homogene sur : Gamma6
%  - Condition de Neumann non homogene sur : Gamma1, Gamma2, Gamma4, Gamma5
%  - Condition de Dirichlet homogene sur: Gamma3
FF  = zeros(NbNodes,1);    % vecteur source du probleme
GG1  = zeros(NbNodes,1);   % vecteur cond. neumann Gamma1
GG2  = zeros(NbNodes,1);   % vecteur cond. neumann Gamma2
GG4  = zeros(NbNodes,1);   % vecteur cond. neumann Gamma4
GG5  = zeros(NbNodes,1);   % vecteur cond. neumann Gamma5
% assemblages
for i=1:NbNodes
    xi=AllNodes(i,1);
    yi=AllNodes(i,2);
    FF(i,1) = FE_f(xi,yi);
    % Gamma1 ==> Neumann non-homogene flux entrant
    if RefNodes(i)==1
        GG1(i,1) = 5*yi; % lineaire en yi: vaut 1 si yi=1/5, vaut 0 si yi=0;
        GG2(i,1) = 0;
        GG4(i,1) = 0;
        GG5(i,1) = 0;
    % Gamma2 ==> Neumann non-homogene flux entrant
    elseif RefNodes(i)==2
        GG1(i,1) = 0;
        GG2(i,1) = -4*xi + 1; % lineare en xi: vaut 1 si xi=0, vaut 0 si xi=1/4
        GG4(i,1) = 0;
        GG5(i,1) = 0;
    % GAmma4 ==> Neumann non-homogene flux sortant
    elseif RefNodes(i)==4
        GG1(i,1) = 0;
        GG2(i,1) = 0;
        GG4(i,1) = -4*xi + 3; %lineare en xi: vaut -1 si xi=1, vaut 0 si xi=3/4
        GG5(i,1) = 0;
    % GAmma5 ==> Neumann non-homogene flux sortant
    elseif RefNodes(i)==5
        GG1(i,1) = 0;
        GG2(i,1) = 0;
        GG4(i,1) = 0;
        GG5(i,1) = -5*yi ; %lineare en yi: vaut -1 si yi=1/5, vaut 0 si yi=0
    else % Neumann homogene (pas de flux)
        GG1(i,1) = 0;
        GG2(i,1) = 0;
        GG4(i,1) = 0;
        GG5(i,1) = 0;
    end
end

LL = zeros(NbNodes,1);  % vecteur second membre
LL0 = zeros(NbNodes,1);% L0
LL1 = zeros(NbNodes,1);% L1
LL2 = zeros(NbNodes,1);% L2
PHIr = FE_relevement(mesh); % relevement
LL = MM*FF + SS1*GG1 + SS2*GG2 + SS4*GG4 + SS5*GG5 - AA*PHIr;
LL0 = MM*FF + SS1*GG1 + SS2*GG2 + SS4*GG4 + SS5*GG5 - AA0*PHIr;
LL1 = - AA1*PHIr;
LL2 = - AA2*PHIr;


% Conditions de Dirichlet homogenes sur Gamma3
% -------------------------------------------
% Note: les noeuds de Dirichlet sont les noeuds du bord Gamma3
Dirichlet_Nodes= find(RefNodes==3);
NbDirichlet_Nodes = size(Dirichlet_Nodes,1);
disp(sprintf('Nombre de noeuds de Dirichlet = %i', NbDirichlet_Nodes));
% Degres de libertes = ensemble des noeuds, a l'exception des 
%                      noeuds de dirichlet
DofNodes = setdiff(transpose(1:NbNodes),Dirichlet_Nodes);
disp(sprintf('Nombre degres de liberte = %i', size(DofNodes,1)));

% extraction des lignes/colonnes qui correspondent aux degres de liberte
AA = AA(DofNodes,DofNodes);
BB = BB(DofNodes,DofNodes);
AA0 = AA0(DofNodes,DofNodes);
AA1 = AA1(DofNodes,DofNodes);
AA2 = AA2(DofNodes,DofNodes);
% second membre : 
% extraction des lignes qui correspondent aux degres de liberte
LL = LL(DofNodes);
LL0 = LL0(DofNodes);
LL1 = LL1(DofNodes);
LL2 = LL2(DofNodes);

AA_decomp{1} = AA0;
AA_decomp{2} = AA1;
AA_decomp{3} = AA2;
LL_decomp{1} = LL0;
LL_decomp{2} = LL1;
LL_decomp{3} = LL2;
end

