function [ gdir ] = FE_gdir( x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE_gdir:
%          
% INPUT * x,y : coordonnees d'un point
%
% OUTPUT * gdir = valeur de gdir(x,y) cond. dirichlet non-homogene
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gdir= 0.6 - 1.2*x;

end

