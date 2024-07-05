function [ kappa ] = FE_kappa_coeff( ref )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUT * ref = valeur de l'indice du sous-domaine (0,1 ou 2)
%
% OUTPUT - kappa= valeur du coefficient de diffusion kappa dans le
%                 sous-domaine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch (ref)
    case 0
        kappa = 1; % valeur fixee
    case 1
        kappa = 0.1; % parametre 1
    case 2
        kappa = 1; % parametre 2
    otherwise
        error('la reference doit etre 0,1 ou 2 !!!')
end
end

