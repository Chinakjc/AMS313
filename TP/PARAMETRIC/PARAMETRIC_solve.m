function [ UU ] = PARAMETRIC_solve( mu, AA_decomp, LL_decomp )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRIC_solve
%          
% INPUT - mu (taille 2) : vecteur des parametres
%       - AA_decomp (cell Qa x 1) : decomposition affine de la matrice AA
%       - LL_decomp (cell Ql x 1) : decomposition affine du second membre LL
%
% OUTPUT - UU (taille NbDof) : solution Haute Fidelite
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nombre de termes decomp affines 
Qa = length(AA_decomp);
Ql = length(LL_decomp);

% taille de l'espace d'approximation Haute fidelite
NbDof = size(AA_decomp{1},1);


AA = sparse(NbDof, NbDof);
LL = zeros(NbDof,1); 
thetaA = PARAMETRIC_thetaA(mu);
thetaL = PARAMETRIC_thetaL(mu);
for q=1:3
        AA = AA + thetaA(q)*AA_decomp{q};
        LL = LL + thetaL(q)*LL_decomp{q};
end

UU = AA \ LL;
   
%error('PARAMETRIC_solve() not yet implemented')


end

