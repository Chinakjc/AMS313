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

avec_AA_decomp = true; 
avec_LL_decomp = true;

% -------------------------
% construction du maillage
% -------------------------
nx = 200;
ny = 40;
mesh = MESH_build_cartesian(nx, ny);

% ---------------------------
% assemblages des matrices EF
% ---------------------------
if (avec_AA_decomp && ~avec_LL_decomp)
    [ DofNodes, AA_ref, LL_ref, ... 
      MM, DDX, DDY, BB, AA_decomp] = FE_assemblages(mesh);
elseif((avec_AA_decomp && avec_LL_decomp)||...
       (~avec_AA_decomp && avec_LL_decomp))
    [ DofNodes, AA_ref, LL_ref,...
      MM, DDX, DDY, BB, AA_decomp, LL_decomp] = FE_assemblages(mesh);
else
    error('il ny a rien a valider: revoir les options utilisateur !')
end

disp('--------------------------------------------');
disp(' Assemblage avec les decompositions affines ');
disp('--------------------------------------------');

% assemblage de la matrice AA
% ---------------------------
NbDof = size(DofNodes,1);
AA = sparse(NbDof, NbDof); 
if (avec_AA_decomp)
    % coeff dans la decomposition affine 
    kappa1 = FE_kappa_coeff(1);
    kappa2 = FE_kappa_coeff(2);
    mu = [kappa1; kappa2];
    thetaA = PARAMETRIC_thetaA(mu);
    % Assemblage a l'aide de la decomposition affine
    disp(' => Assemblage matrice')
    for q=1:3
        AA = AA + thetaA(q)*AA_decomp{q};
    end
else
    % sinon on recupere la reference
    AA = AA_ref;
end


% assemblage du second membre LL
% ---------------------------
LL = zeros(NbDof,1); 
if (avec_LL_decomp)
    % coeff dans la decomposition affine 
    kappa1 = FE_kappa_coeff(1);
    kappa2 = FE_kappa_coeff(2);
    mu = [kappa1; kappa2];
    thetaL = PARAMETRIC_thetaL(mu);
    % Assemblage a l'aide de la decomposition affine
    disp(' => Assemblage second membre')
    for q=1:3
        LL = LL + thetaL(q)*LL_decomp{q};
    end
else
    % sinon on recupere la reference
    LL = LL_ref;
end


disp('-------------------------');
disp(' Resolution des systemes ');
disp('-------------------------');
% solution de reference (= sans utiliser les decompositions affines)
disp(' => syteme de reference')
UU_ref = AA_ref\ LL_ref;
disp(' => syteme avec decomposition affine')
UU = AA \ LL;

% calcul de l'erreur
% ------------------
% Note: lorsque tout est correct, on doit trouver la precision machine
diff = UU - UU_ref;
disp(sprintf('Erreur absolue = %e',norm(diff)))
disp(sprintf('Erreur relative= %e',norm(diff)/ norm(UU_ref)))