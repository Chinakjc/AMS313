function  UU_perp  = LINALG_orthonormalize( UU, PP, BB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINALG_orthonormalize :
% orthonormaliser le vecteur UU par rapport au produit scalaire BB
%          
% INPUT * UU: vecteur NbNodes a orthogonaliser dans la base PP
%       * PP (taille NbNodes x N) : matrice des N vecteurs de base contre
%       lesquels il faut orthonormaliser
%       * BB: matrice du produit scalaire
%
% OUTPUT - UU_perp: vecteur orthonormalise,
%                   UU_perp = UU - Proj_PP[UU]
%                   ou Proj_PP[] est le projecteur orthonormal sur PP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cas ou P est vide
if length(PP)== 0
    normUU = sqrt(UU'*BB*UU);
    UU_perp = UU / normUU;
% cas general
else
    % nombre de vecteur contre lesquels orthogonaliser
    N = size(PP,2);
    
    normUU = sqrt(UU'*BB*UU);
    UU_perp = UU / normUU;
    
    normUU_perp= 0.0;
    
    %gram-schmidt with re-iteration
    max_iter= 10;
    iter=0;
    while (normUU_perp < 0.1)&&(iter < max_iter)
        % orthogonaliser en
        % soustrayant la projection sur les vecteurs de bases existants
        for i=1:N
            UU_perp = UU_perp - (PP(:,i)'*BB*UU_perp)*PP(:,i) ; 
        end

        % normaliser
        normUU_perp = sqrt(UU_perp'*BB*UU_perp);
        UU_perp = UU_perp / normUU_perp;
        
        % maj
        iter=iter+1;
    end
    if(iter==max_iter)
        disp('Gram-Schmidt with re-iteration: found colinear vectors!')
    end
end

end

