% ============================================================
% Construction d'une base reduite par la methode POD 
% pour l'equation de Poisson 2D, avec conditions de
% Neumann et de Dirichlet sur le bord
%  
% avec 2 parametres:
%   * kappa1 (coeff diffusion dans ss-domaine \Omega1)
%   * kappa2 (coeff diffusion dans ss-domaine \Omega2)
%
% Note: le coeff diffusion dans le ss-domaine \Omega0 est fixe
% ============================================================

close all;
clear all;

% -------------------------
% definitions 
% -------------------------
% bornes du domaine parametrique
kappa1_min = 0.1;
kappa1_max = 1;
kappa2_min = 0.1;
kappa2_max = 1;
type_plan_experience = 'cartesien';%'random';%'cartesien';%'random';%'cartesien';%'cartesien'; % 'cartesien' ou 'random'
sauvegarde_base_reduite = true;    % true pour sauvegarder la base reduite
                                    % PP dans un fichier "PP.mat"

% -------------------------
% construction du maillage
% -------------------------
nx = 200;
ny = 40;
mesh = MESH_build_cartesian(nx, ny);

% assemblages des matrices EF
% ---------------------------
[ DofNodes, AA_ref, LL_ref,...
      MM, DDX, DDY, BB, AA_decomp, LL_decomp] = FE_assemblages(mesh);
  
% nombre de degres de liberte
NbDof = size(DofNodes,1);

disp('--------------------');
disp(' Phase Exploratoire ');
disp('--------------------');
tic;

% definition du plan d'experience 
% --------------------------------
% on construit tous les couples (kappa1, kappa2)
% pour lesquels le probleme EF sera resolu (n_train couples en tout)
mu_list = []; % (taille n_train x 2)

if (strcmp(type_plan_experience, 'cartesien'))
    size_train = 12;
    mu_list_1 = linspace(kappa1_min,kappa1_max,size_train);
    mu_list_2 = linspace(kappa2_min,kappa2_max,size_train);
    [MU1,MU2] = meshgrid(mu_list_1,mu_list_2);
    MU1V = reshape(MU1,[],1);
    MU2V = reshape(MU2,[],1);
    mu_list =  [MU1V,MU2V];
    %error('type_plan_experience cartesien not yet implemented')
elseif (strcmp(type_plan_experience,'random'))
    size_train = 35;
    scale1 = kappa1_max - kappa1_min;
    scale2 = kappa2_max - kappa2_min;
    scale = [scale1, 0; 0, scale2];
    offset = [kappa1_min, kappa2_min];
    for  q = 1:size_train
        mu_list = [mu_list ; offset + rand(1,2) * scale];
    end
    %error('type_plan_experience random not yet implemented')
else
    error('type_plan_experience pas bien defini')
end

% resolution EF pour tout les couples (kappa1, kappa2) 
% -------------------------------------------------
n_train = size(mu_list,1);
disp(sprintf('Nbre de resolutions HF = %i', n_train));
AllUUs = zeros(NbDof,n_train);
for j=1:n_train
     % valeur du parametre
    mu_j = mu_list(j,:); 
    % calcul de la solution HF
    UU = PARAMETRIC_solve( mu_j, AA_decomp, LL_decomp );
    % on sauvegarde la solution dans la colonne j de AllUUs
    AllUUs(:,j) = UU;
end

elapsed = toc;
disp(sprintf('Phase Exploration elapsed time = %f s', elapsed));

disp('-------------------');
disp(' Phase Compression ');
disp('-------------------');
tic;

CVM = AllUUs' * BB   / n_train * AllUUs;
elapsed = toc;
disp(sprintf('Phase Compression elapsed time = %f s', elapsed));

% Courbes de decroissance des valeurs propres
% ------------------------------------------------
[V,D] = eig(CVM);
eigenvalues = real(diag(D));
ZZ = AllUUs * real(V);

[sortedEigenvalues, sortOrder] = sort(eigenvalues, 'descend');


plot(1:size(sortedEigenvalues),sortedEigenvalues, '-o');  
xlabel('Index');  
ylabel('Eigenvalue');  
set(gca, 'YScale', 'log');
title('Eigenvalues from Largest to Smallest');  

% Troncature et definition d'une base reduite
% -------------------------------------------
affichage = false;
for N = [4,8,12,16,20,24,28,32]
    PP = []; % matrice (taille NbDof x N) representant une base reduite de taille N
    for q = 1:N%sortOrder(1:N)
        PP = [PP, real(LINALG_orthonormalize( ZZ(:,q), PP, BB))];
    end

%size(PP)

    if(affichage)
        for q = 1:min(N,5)
            FE_visu(FE_add_Dirichlet_DoFs(PP(:,q), mesh , DofNodes), mesh, sprintf('fonction de base reduite %d',q));
        end
    end


% Sauvergarde de la base reduite construite
% -----------------------------------------
    if(sauvegarde_base_reduite)
        disp('Sauvegarde de la base reduite dans le fichier PP_pod.mat');
        save(sprintf("PP_pod_%d.mat",N), "PP");
    end

end