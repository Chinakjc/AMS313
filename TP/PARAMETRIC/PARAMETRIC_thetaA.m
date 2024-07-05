function [ thetaA ] = PARAMETRIC_thetaA( mu )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRIC_thetaA :
% coefficient dans la decomposition affine de AA
%          
% INPUT * mu (taille 2) parametre

% OUTPUT - thetaA (taille 3) vecteur des coefficients de la decomposition
%           affine de AA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa1 = mu(1);
kappa2 = mu(2);
thetaA = [1; kappa1; kappa2];

end

