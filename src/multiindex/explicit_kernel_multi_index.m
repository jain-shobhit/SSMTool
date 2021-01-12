function [K_k,G_k,innerresonance] = explicit_kernel_multi_index(z_k, Lambda_M_vector, Lambda_Mk_vector, W_M, reltol)
%% COEFF_MATR_KERNEL Explicit kernel-construction of the coefficient-matrix
% This function computes the kernel of the coefficient matrix for eigenvalue 
% pairs that are in resonance as described in the document ''Explicit Kernel Extraction 
% and Proof ofSymmetries of SSM Coefficients - Multi-Indexversion''.
%SSM dimension
l         = size(Lambda_M_vector,1);
%Compare for all combinations if singularity occurs
Lambda_Ci = Lambda_M_vector - Lambda_Mk_vector; % column vector - row vector
%threshold below which resonance occurs
ref       = min(abs(Lambda_M_vector));
abstol = reltol*ref;
%find eigenvalues that trigger resonance
[I,F]  = find(abs(Lambda_Ci) < abstol); % I for eigenvalue and F for combination
r_k = length(I);
if r_k
    innerresonance = 1;
    
    % create E_F, E_I
    E_F = sparse( F, (1:r_k).', true(r_k,1), z_k, r_k);
    E_I = sparse( I, (1:r_k).', true(r_k,1), l, r_k);
    
    % create K_k, G_k
    K_k = khatri_rao_product(E_F, W_M(:,I));
    G_k = khatri_rao_product(E_F, E_I)';
else
    innerresonance = 0;
    
    K_k=[];
    G_k=[];
end