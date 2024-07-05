function [X] = FE_solve(mu, AA_decomp, LL_decomp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
NbDof = size(LL_decomp{1},1);
AA = sparse(NbDof, NbDof); 
LL = zeros(NbDof,1); 
thetaA = PARAMETRIC_thetaA(mu);
thetaL = PARAMETRIC_thetaL(mu);


for q=1:3
    AA = AA + thetaA(q)*AA_decomp{q};
    LL = LL + thetaL(q)*LL_decomp{q};
end


X = AA \ LL; 