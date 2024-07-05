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

max_errors = [];
mean_errors = [];
max_errors2 = [];

% test ramdon
% ---------------------------
size_train = 5;
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

Ns = [4,8,12,16,20,24,28,32];
for N = Ns
    load(sprintf("PP_pod_%d.mat",N))
    [Arb_decomp, Lrb_decomp] = RB_reduced_decomp(AA_decomp,LL_decomp, PP);

    %size(Lrb_decomp{1})



    errors = [];
    errors2 = [];
    for mu = mu_list'  
        Xrb = RB_solve(mu, Arb_decomp, Lrb_decomp);
        Xref = FE_solve(mu, AA_decomp, LL_decomp);
        %FE_visu(FE_add_Dirichlet_DoFs(PP * Xrb, mesh , DofNodes), ...
        %mesh, sprintf('sol RB with mu1 = %f, mu2 = %f',mu(1),mu(2)));
        dX = PP * Xrb - Xref;
        pU = PP * PP' * BB * Xref;
        dU = Xref - pU;
        error = dX' * BB * dX;
        errors = [errors, error];
        error2 = dU' * BB * dU;
        errors2 = [errors2,error2];
    end
    max_error = max(errors);
    mean_error = sum(errors)/size(mu_list,1);
    max_error2 = max(errors2);

    max_errors = [max_errors,max_error];
    mean_errors = [mean_errors,mean_error];
    max_errors2 = [max_errors2,max_error2];

    fprintf("max error = %f\n",max_error);
    fprintf("mean error = %f\n",mean_error);
    fprintf("max error 2 = %f\n",max_error2);

    clear PP;

end


figure(1);
hold on;
plot(Ns,log10(max_errors),'-o');
plot(Ns,log10(mean_errors),'-o');
plot(Ns,log10(max_errors2),'-o');
grid on;
legend('max error (log10)','mean error (log10)','max error 2 (log10)');


