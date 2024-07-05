function [ thetaL ] = PARAMETRIC_thetaL( mu )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRIC_thetaL :
% coefficient dans la decomposition affine de LL
%          
% INPUT * mu (taille 2) parametre

% OUTPUT - thetaA (taille 3) vecteur des coefficients de la decomposition
%           affine de LL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa1 = mu(1);
kappa2 = mu(2);
thetaL = [1; kappa1; kappa2];

end

