% ============================================================
% Version parametrique du solveur elements finis 
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

% -----------------------------------
% Options a definir par l'utilisateur
% -----------------------------------
% Note: les valeurs des parametres kappa1, kappa2 utilisees pour la validation 
% sont celles definies dans FE_kappa_coeff



% -------------------------
% construction du maillage
% -------------------------
nx = 200;
ny = 40;
mesh = MESH_build_cartesian(nx, ny);

% ---------------------------
% assemblages des matrices EF
% ---------------------------

[ DofNodes, AA_ref, LL_ref, ... 
      MM, DDX, DDY, BB, AA_decomp,LL_decomp] = FE_assemblages(mesh);


disp('--------------------------------------------');
disp(' Assemblage avec les decompositions affines RB');
disp('--------------------------------------------');

% assemblage 
% ---------------------------
NbDof = size(DofNodes,1);



% test ramdon
% ---------------------------
size_train = 64;
kappa1_min = 0.1;
kappa1_max = 1;
kappa2_min = 0.1;
kappa2_max = 1;
scale1 = kappa1_max - kappa1_min;
scale2 = kappa2_max - kappa2_min;
scale = [scale1, 0; 0, scale2];
offset = [kappa1_min, kappa2_min];
mu_list = [];

for  q = 1:size_train
    mu_list = [mu_list ; offset + rand(1,2) * scale];
end


load("PP_greedy.mat")
[Arb_decomp, Lrb_decomp] = RB_reduced_decomp(AA_decomp,LL_decomp, PP);
load("QRF.mat")

errors = [];
residus = [];
for mu = mu_list'  
    Xrb = RB_solve(mu, Arb_decomp, Lrb_decomp);
    Xref = FE_solve(mu, AA_decomp, LL_decomp);
        
    dX = PP * Xrb - Xref;
    pU = PP * PP' * BB * Xref;
    dU = Xref - pU;
    error2 = dX' * BB * dX;
    errors = [errors, sqrt(error2)];
    residu2 = RB_compute_resnorm2_stabilized(mu,Xrb,Rmat,f_parallel,f_perp);
    residus = [residus,sqrt(residu2)];
        
end
max_error = max(errors);
mean_error = sum(errors)/size(mu_list,1);
    

  

fprintf("log10 max error = %f\n",log10(max_error));
fprintf("log10 mean error = %f\n",log10(mean_error));
fprintf("log10 max residu = %f\n",log10(max(residus)));
fprintf("log10 mean residu = %f\n",log10(sum(residus)/size(mu_list,1)));

figure
scatter(mu_list(:,1), mu_list(:,2), 100, log10(residus./errors), 'filled'); 
colorbar; 
title("log10 (residus/errors) ");
figure
scatter(mu_list(:,1), mu_list(:,2), 100, log10(residus), 'filled'); 
colorbar; 
title("log10 (residus) ");
figure
scatter(mu_list(:,1), mu_list(:,2), 100, log10(errors), 'filled'); 
colorbar; 
title("log10 (errors) ");

  
    

clear PP;
