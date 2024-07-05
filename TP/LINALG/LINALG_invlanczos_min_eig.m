function [ lambda_min, v , conv_curve, it] = LINALG_invlanczos_min_eig( AA, XX, u1, tol )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Smallest eigenvalue \lambda such that:
%
%         AA v = \lambda XX v,
%
%   with normalisation :
%         v^{T} XX v = 1
%
%   Input: * mu: scalar value of parameter
%          * AA: sparse N x N matrix 
%          * XX: sparse N x N matrix (hermitian, positive definite)
%          * u1: starting vector size N (optional)
%          * tol: desired accuracy (optional)
%
%   Output: * lambda_min: smallest eigenvalue
%           * v: associated eigenvector
%           * conv_curve: error indicator throughout convergence
%           * it: number of iterations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nbpt = size(XX, 1);

% case where u1 and tol are default
if nargin ==2
    % define a (normalized) default starting vector
    u1 = rand(Nbpt,1);
    % normalize:
    % Solve linear system: X (invX_u1) = u1
    invX_u1 = XX \ u1;
    normu1_2 = u1'*invX_u1;
    u1 = u1 / sqrt(normu1_2);
    % define default tolerance
    tol= 1e-9;       
end

% case where tol is default
if nargin ==3
    tol= 1e-9;       
end

% initialize loop
curr_err = 1e9;
conv_curve = [];
maxit = 100;
it = 0;
% storage tridiagonal matrix T
alphas = [];
betas = [];
%
u = u1;
p = 0;
Qmat = [];
% Solve X r = u for r
r = XX \ u;
beta = sqrt(u'*r);
%
while (curr_err > tol)&&(it < maxit)
    %
    q = u / beta;
    Qmat = [Qmat q];
    u_bar = AA\ q;
    u_bar = u_bar - p*beta;
    alpha = q'*u_bar;
    % add to tridiagonal matrix
    alphas= [alphas; alpha];
    %
    p = r / beta;
    r = u_bar - p*alpha;
    u = XX*r;
    %
    % eigenvalue problem here
    [Smat, thetas] = LINALG_tridiagonal_eigensolve(alphas, betas);
    %
    beta = sqrt(u'*r);
    %
    % error estimation (on the largest eigenvalue):
    % ------
    % find index of largest Ritz value
    [theta_max, ind] = max(thetas);
    curr_err = beta*abs(Smat(it+1, ind));
    conv_curve = [conv_curve; curr_err];
    % ------
    %
    % add to tridiagonal matrix
    betas= [betas; beta];
    %
    it = it + 1;
end

% the smallest eigenvalue is well approximated by the inverse of the largest Ritz value:
lambda_min = 1 / theta_max;
% the associated eigenvector is the Ritz vector
vhat = Qmat*Smat(:, ind);
v = XX \ vhat;

end



