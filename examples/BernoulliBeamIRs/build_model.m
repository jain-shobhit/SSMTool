function [M,C,K,fnl] = build_model(nElements,kLinear,kNonlinear)

[M,C,K]=L_Bernoulli_Beam_Model(nElements);
n = length(M);
K(n-1,n-1) = K(n-1,n-1)+kLinear; 

%% 
% first-order tensors

f2 = sptensor([n,n,n]);
f3 = sptensor([n,n,n,n]);
%% 
% adding cubic spring to the end node of the beam
dof = n-1;
for j = dof
    f3(j,j,j,j) = kNonlinear;
end

fnl = {f2,f3};