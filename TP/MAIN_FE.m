% =====================================================
% Solveur elements finis 
% pour l'equation de Poisson 2D suivante, avec conditions de
% Neumann et de Dirichlet sur le bord
%
% | -div (K grad \phi(x,y)) = f(x,y),          dans \Omega
% |               \phi(x,y) = gdir(x,y),       sur G3 
% |      K grad \phi(x,y).n = gneu(x,y),       sur G1 U G2 U G4 U G5,
% |        grad \phi(x,y).n = 0,               sur G6
% 
%      G2    G3    G4
%    |----|------|----|
% G1 |     \Omega     | G5
%    |----------------|
%            G6
% 
% K = coeff de diffusion (cf. FE_kappa_coeff)
% gdir = donnee de dirichlet non-homogene (cf. FE_gdir)
% gneu = donnee de neumann non-homogene (cf. FE_assemblages)
% f = second membre (cf. FE_f)
% =====================================================

close all;
clear all;

% -----------------------------------
% Options a definir par l'utilisateur
% -----------------------------------
calcul_gradient_phi = true;

% -------------------------
% construction du maillage
% -------------------------
nx = 200;
ny = 40;
mesh = MESH_build_cartesian(nx, ny);

% assemblages des matrices EF
% ---------------------------
[ DofNodes, AA, LL, MM, DDX, DDY] = FE_assemblages(mesh);
%[ DofNodes, AA_ref, LL_ref,...
%      MM, DDX, DDY, BB, AA_decomp, LL_decomp] = FE_assemblages(mesh);

disp('-----------------------');
disp(' Resolution du systeme ');
disp('-----------------------');
UU = AA \LL; 
%mu = [FE_kappa_coeff(1); FE_kappa_coeff(2)];
%UU = PARAMETRIC_solve( mu, AA_decomp, LL_decomp );

% ajout des noeuds de dirichlet
UU_full = FE_add_Dirichlet_DoFs( UU, mesh , DofNodes);

% prise en compte du relevement
PHIr = FE_relevement(mesh);
PHI = UU_full + PHIr;

% visualisation de la solution 
% -----------------------------
FE_visu(UU_full, mesh, 'Solution U');
FE_visu(PHI, mesh, 'Solution \Phi');

% Approx des derivees de PHI /x et /y
% -----------------------------------
if (calcul_gradient_phi)
    % approximation de la derivee de la solution /x
    LLX = -DDX*PHI;
    PHI_DX = MM \ LLX;
    % approximation de la derivee de la solution /y
    LLY = -DDY*PHI;
    PHI_DY = MM \ LLY;
    % affichages
    FE_visu(PHI_DX, mesh, 'Derivee d\Phi / dx ');
    FE_visu(PHI_DY, mesh, 'Derivee d\Phi / dy ');
    FE_visu_vector(PHI_DX, PHI_DY,mesh, 'Gradient de \Phi')
end