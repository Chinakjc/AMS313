function [Smat, thetas] = NNSCM_tridiagonal_eigensolve(alphas, betas)
%   Eigenvalues and eigenvectors of a tridiagonal matrix
%
%   Input: * alphas: diagonal terms (size m)
%          * betas: extra-diagonal terms (size m-1)
%
%   Output: * Smat: matrix of size m x m containg eigenvectors
%           * thetas: eigenvalues
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = length(alphas);
m_minus_one = length(betas);
if m_minus_one ~= m - 1
    disp('ERROR: alphas and betas not coherent to form a tridiagonal matrix')
    exit
end

% build tridiagonal matrix;
if m==1
    % special case where betas is empty
    T = diag(alphas);
else
    T = diag(alphas) + diag(betas, 1) + diag(betas, -1);
end

% solve eigenvalue problem
[Smat, D] = eig(T);

thetas = diag(D);

end

