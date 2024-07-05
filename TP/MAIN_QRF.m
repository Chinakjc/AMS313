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

load("PP_greedy_app.mat")

% construction du maillage
% -------------------------
nx = 200;
ny = 40;
mesh = MESH_build_cartesian(nx, ny);

% assemblages des matrices EF
% ---------------------------
[ DofNodes, AA_ref, LL_ref,...
      MM, DDX, DDY, BB, AA_decomp, LL_decomp] = FE_assemblages(mesh);


N = size(PP,2);
% n=1,...,N (taille Nbpt x N*Qa)
Qmat = [];
% matrice de changement de base (taille N*Qa x N*Qa)
Rmat = [];
% partie du RHS qui appartient au sous espace Range(Qmat), taille Qa*N x Ql
f_parallel=[]; 

Qa = 3;
Ql = 3;

LL_hat_decomp = cell(Ql,1);
for q=1:Ql
    LL_hat_decomp{q} = BB \ LL_decomp{q};
end

for q = 1:N
    PP_q = PP(:,1:q);
    new_basis_func = PP(:,q);
    [Arb_decomp, Lrb_decomp] = RB_reduced_decomp(AA_decomp,LL_decomp, PP);
    [Qmat, Rmat]= RB_update_Qmat_Rmat(N, new_basis_func, BB, AA_decomp, Qmat, Rmat);
    [f_parallel, f_perp] = RB_update_resparts_stabilized(Qa, BB, Qmat, LL_decomp, LL_hat_decomp,f_parallel);
end

disp('Sauvegarde de Qmat, Rmat, f_parallel, f_perp');
save('QRF.mat', 'Qmat','Rmat','f_parallel','f_perp');