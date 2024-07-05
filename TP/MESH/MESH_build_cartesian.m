function [ mesh ] = MESH_build_cartesian( nx, ny)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction maillage cartesien pour le domaine \Omega suivant:
%
%      G2    G3    G4
%    |----|------|----|
% G1 |     \Omega     | G5
%    |----------------|
%            G6
%
% INPUTS : nx (integer) = nombre de maille selon x 
%          ny (integer) = nombre de maille selon y
%
% OUTPUT : mesh (structure) = tous les tableaux concernant le maillage
%
% La structure mesh contient les attributs (liste non exhaustive) :
%       * NbNodes : nombre de noeuds
%       * NbTriangles : nombre de triangles
%       * NbBoundaryEdges : nombre d'aretes de bord
%       * AllNodes (taille NbNodes x 2): coordonnees x,y des noeuds
%       * AllTriangles (taille NbTriangles x 3): indices 3 noeuds de chaque 
%         triangle
%       * RefNodes (taille NbNodes) : references des noeuds
%       * RefTriangles (taille NbTriangles) : references des triangles
%       * AllBoundaryEdges (taille NbBoundaryEdges x 2): indices 2 noeuds
%         de chaque arete de bord
%       * RefBoundaryEdges (taille NbBoundaryEdges): references des aretes
%         de bord
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialisation topologie
MESH_global_params;

% reference des noeuds/aretes du bord
ref_bas = 6;
ref_gauche = 1;
ref_droite = 5;
ref_haut = 2;

% taille du domaine
xmin = 0.;
xmax = 1;
ymin = 0.;
ymax = 0.2;

% Definition des mailles quadrangles 
%                     et triangles
% ----------------------------------
% 1 quad = 4 indices de noeud
% 1 triangle = 3 indices de noeud
AllQuads = zeros(nx*ny, 4);
AllTriangles = zeros(2*nx*ny, 3);
for j=1:ny
    for i=1:nx
        % indice du quad
        q = (j-1)*nx + i;
        % indice du noeud en bas a gauche
        AllQuads(q,topology.bottomleft_ind) = q + j - 1;
        % indice du noeud en bas a droite
        AllQuads(q,topology.bottomright_ind) = q + j;
        % indice du noeud en haut a gauche
        AllQuads(q,topology.topleft_ind) = (j*nx + i) + j;
        % indice du noeud en haut a droite
        AllQuads(q,topology.topright_ind) = (j*nx + i) + j + 1;
        % construction de 2 triangles a partir du quadrangle
        t1 = 2*(q-1) +1;
        t2 = 2*q;
        % decoupe du quadrangle selon l'une ou l'autre diagonale, au hasard
        if (rand(1,1)> 0.5)
            % decoupage selon la diagonale qui va d'en bas a gauche vers 
            % en haut a droite
            AllTriangles(t1,1) = AllQuads(q,1);
            AllTriangles(t1,2) = AllQuads(q,2);
            AllTriangles(t1,3) = AllQuads(q,4); 
            AllTriangles(t2,1) = AllQuads(q,1);
            AllTriangles(t2,2) = AllQuads(q,3);
            AllTriangles(t2,3) = AllQuads(q,4); 
        else
            % decoupage selon la diagonale qui va d'en haut a gauche vers 
            % en bas a droite
            AllTriangles(t1,1) = AllQuads(q,1);
            AllTriangles(t1,2) = AllQuads(q,2);
            AllTriangles(t1,3) = AllQuads(q,3);
            AllTriangles(t2,1) = AllQuads(q,2);
            AllTriangles(t2,2) = AllQuads(q,3);
            AllTriangles(t2,3) = AllQuads(q,4); 
        end
    end
end

% Definition des voisins de chaque quad
% ------------------------------------- 
% NeighQuads(q,:)= 4 indices des 4 quad voisins du quad q
% N.B. -1 si pas de voisin (quad au bord du domaine)
NeighQuads = -ones(nx*ny, 4);
for j=1:ny
    for i=1:nx
        % indice du quad
        q = (j-1)*nx + i;
        % voisin de gauche (pas de voisin de gauche si i==1)
        if (i~=1) 
            NeighQuads(q,topology.left_ind) = q-1;
        end
        % voisin du bas (pas de voisin du bas si j==1)
        if (j~=1) 
            NeighQuads(q,topology.bottom_ind) = q-nx;
        end
        % voisin de droite (pas de voisin de droite si i==nx)
        if (i~=nx) 
            NeighQuads(q,topology.right_ind) = q+1;
        end
        % voisin du haut (pas de voisin de droite si j==ny)
        if (j~=ny) 
            NeighQuads(q,topology.top_ind) = q+nx;
        end
    end
end

% Definition des noeuds 
% -----------------------
% 1 noeud = 2 coordonnees x, y
Nbnodes = (nx + 1)*(ny + 1);
AllNodes = zeros(Nbnodes, 2);
RefNodes = zeros(Nbnodes, 1);
% taille des mailles selon x et selon y
hx = (xmax-xmin)/nx;
hy = (ymax-ymin)/ny;
for j=1:(ny+1)
    yj = ymin + (j-1)*hy;
    for i=1:(nx+1)
        xi = xmin + (i-1)*hx;
        n = (j-1)*(nx+1) + i;
        AllNodes(n,1) = xi;
        AllNodes(n,2) = yj;
        % Bord du Bas
        if (j==1)
            RefNodes(n) = ref_bas;
        end
        % Bord du Haut 
        if (j==ny+1)
            RefNodes(n) = ref_haut;
        end
        % Bord de Gauche
        if (i==1)
            RefNodes(n) = ref_gauche;
        end
        % Bord de Droite
        if (i==nx+1)
            RefNodes(n) = ref_droite;
        end
    end
end

% Definition des aretes de bord
% ----------------------------------
% 1 arete = 2 indices de noeud
Nbedges = 2*(nx+ny);
AllBoundaryEdges = zeros(Nbedges, 2);
RefBoundaryEdges = zeros(Nbedges, 1);
% bords haut et bas
for i=1:nx
    % bord bas
    % indice
    a = 2*(i-1)+1;
    RefBoundaryEdges(a) = ref_bas;
    q = i;
    AllBoundaryEdges(a,1)= AllQuads(q,topology.bottomleft_ind); % indice noeud en bas a gauche
    AllBoundaryEdges(a,2)= AllQuads(q,topology.bottomright_ind); % indice noeud en bas a droite
    % bord haut
    % indice
    a = 2*i;
    RefBoundaryEdges(a) = ref_haut;
    q = (ny-1)*nx + i;
    AllBoundaryEdges(a,1)= AllQuads(q,topology.topleft_ind); % indice noeud en haut a gauche
    AllBoundaryEdges(a,2)= AllQuads(q,topology.topright_ind); % indice noeud en haut a droite
end
% bords gauche et droite
for j=1:ny
    % bord gauche
    % indice
    a = 2*nx + 2*(j-1)+1;
    RefBoundaryEdges(a) = ref_gauche;
    q = nx*(j-1) + 1;
    AllBoundaryEdges(a,1) = AllQuads(q,topology.bottomleft_ind); % indice noeud en bas a gauche
    AllBoundaryEdges(a,2) = AllQuads(q,topology.topleft_ind); % indice noeud en haut a gauche
    % bord droit
    % indice
    a = 2*(nx +j);
    RefBoundaryEdges(a) = ref_droite;
    q = nx*j;
    AllBoundaryEdges(a,1) = AllQuads(q,topology.bottomright_ind); % indice noeud en bas a droite
    AllBoundaryEdges(a,2) = AllQuads(q,topology.topright_ind); % indice noeud en haut a droite
end

% -------------------------------
% Cette partie depend du probleme 
% -------------------------------
% 1/ on separe le bord haut en trois bord distincts
% ref=2 entre 0 <= x <= 0.25
% ref=3 entre 0.25 <= x <= 0.75
% ref=4 entre 0.75 <= x <= 1

% traitement des noeuds de bord
nodes_bord_haut = find(RefNodes==ref_haut);
for i=1:size(nodes_bord_haut)
    n = nodes_bord_haut(i);
    x = AllNodes(n,1);
    if (x <= 0.25)
        RefNodes(n) = 2;
    elseif (x <= 0.75)
        RefNodes(n) = 3;
    else
        RefNodes(n) = 4;
    end
end
% traitement des aretes de bord
edges_bord_haut = find(RefBoundaryEdges==ref_haut);
for i=1:size(edges_bord_haut)
    a = edges_bord_haut(i);
    n1 = AllBoundaryEdges(a,1);
    n2 = AllBoundaryEdges(a,2);
    x_moy = 0.5*(AllNodes(n1,1) + AllNodes(n2,1));
    if (x_moy <= 0.25)
        RefBoundaryEdges(a) = 2;
    elseif (x_moy <= 0.75)
        RefBoundaryEdges(a) = 3;
    else
        RefBoundaryEdges(a) = 4;
    end
end

% 2/ Le domaine \Omega est divise en 3 parties
% \Omega1 = ]1/3,4/9[   X ]0,1/5[ (ref==1)
% \Omega2 = ]4/9,5/9[   X ]0,1/5[ (ref==2)
% \Omega3 = ]5/9,2/3[   X ]0,1/5[ (ref==3)
% \Omega0 = \Omega prive de (union \Omegai, i=1,2,3)
RefQuads = zeros(nx*ny, 1);
RefTriangles = zeros(2*nx*ny, 1);
for j=1:ny
    for i=1:nx
        % indice du quad
        q = (j-1)*nx + i;
        % indices des 2 triangles construits a partir du quadrangle
        t1 = 2*(q-1) +1;
        t2 = 2*q;
        % On recupere les coordonnees des 4 noeuds du quad
        S1 = AllNodes(AllQuads(q,topology.bottomright_ind),:);
        S2 = AllNodes(AllQuads(q,topology.topright_ind),:);
        S3 = AllNodes(AllQuads(q,topology.bottomleft_ind),:);
        S4 = AllNodes(AllQuads(q,topology.topleft_ind),:);
        % calcul de l'isobarycentre du quad
        SG = 0.25*(S1 + S2 + S3 + S4);
        % on discrimine selon les coordonnees x,y du barycentre
        x = SG(1);
        y = SG(2);
        ref_q = 0;
        if (4/9 <= x)&&(x <= 5/9)&&(0.06 <= y)&&(y <= 0.14)
            ref_q = 2;
        elseif (1/3 <= x)&&(x <= 2/3)&&(0.02 <= y)&&(y <= 0.18)
            ref_q = 1;
        end
        % on rempli les tableaux avec la bonne reference
        RefQuads(q,1) = ref_q;
        RefTriangles(t1) = ref_q;
        RefTriangles(t2) = ref_q;
    end
end

% On recupere tous les tableaux dans une structure unique
% -------------------------------------------------------
mesh = struct;
mesh.NbTriangles = size(AllTriangles,1);
mesh.NbNodes = size(AllNodes,1);
mesh.NbBoundaryEdges = size(RefBoundaryEdges,1);
mesh.AllQuads = AllQuads;
mesh.AllTriangles = AllTriangles;
mesh.AllNodes = AllNodes;
mesh.RefQuads = RefQuads;
mesh.RefTriangles = RefTriangles;
mesh.RefNodes = RefNodes;
mesh.AllBoundaryEdges = AllBoundaryEdges;
mesh.RefBoundaryEdges = RefBoundaryEdges;
mesh.NeighQuads = NeighQuads;

end

