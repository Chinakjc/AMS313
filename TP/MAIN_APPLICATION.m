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
nx = 400;
ny = 80;
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
JJ = FE_quantity_of_interest( mesh, DofNodes );


% test ramdon
% ---------------------------
size_train = 3000;
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


load("PP_greedy_app.mat")
[Arb_decomp, Lrb_decomp] = RB_reduced_decomp(AA_decomp,LL_decomp, PP);
load("QRF_app.mat")

Jrb = JJ' * PP;

Qrbs = [];

disp('--------------------------------------------');
disp(' Test Starts');
disp('--------------------------------------------');
tStart = tic;

NJ = norm(JJ);

dQs = zeros(size_train,1);
q = 1;
T = zeros(1,size_train);
for mu = mu_list'  
    tic
    Xrb = RB_solve(mu, Arb_decomp, Lrb_decomp);
    Qrb = Jrb * Xrb;
    Qrbs = [Qrbs, Qrb];
    T(q) = toc;
    residu2 = RB_compute_resnorm2_stabilized(mu,Xrb,Rmat,f_parallel,f_perp);
    dQs(q) = NJ * sqrt(residu2);
    q = q + 1;
end
tExp = sum(T)
tEnd = toc(tStart)


figure
scatter(mu_list(:,1), mu_list(:,2), 100, Qrbs, 'filled'); 
colorbar; 
title("Qrb");

figure
scatter(mu_list(:,1), mu_list(:,2), 100, log10(dQs), 'filled'); 
colorbar;
title("log10 (dQ)");