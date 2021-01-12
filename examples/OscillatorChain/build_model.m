function [M,C,K,fnl,f_ext] = build_model(n,m,c,k,kappa2,kappa3)

[K,C,f2,f3] = assemble_global_coefficients(k,kappa2,kappa3,c,n);
%% 
% Dirichlet boundary conditions

K = K(2:n+1,2:n+1);
C = C(2:n+1,2:n+1);
M = m*speye(n,n);
f2 = f2(2:n+1,2:n+1,2:n+1);
f3 = f3(2:n+1,2:n+1,2:n+1,2:n+1);

fnl = {f2, f3};

f_ext = ones(n,1);
% f_ext = 0.1*ones(n,1)+8*(1:n)'/n; % ones(n,1) has no contribution to the second mode
