function[W_0_1,R_0_1, multi_input] = coeffs_setup(obj,order)
%% SSM_MULTI Setup for the mutli-index calculation
% This function checks if the symmetries described in $\textit{Explicit Kernel
% Extraction and Proof ofSymmetries of SSM Coefficients - Multi-Indexversion}$
% can be used. This is the  case if the inputs $\texttt{A,B,F}$ are all purely
% real and the eigenvectors are compl. conjugate for complex conjugate eigenvalues
%, we say the system is real. If this is not the case, then the full coefficients
% get calculated and we say the system is not (purely) real. The calculation is
% executed in $\texttt{SSM\_coefficients}$. This function is designed to distinguish
% between the case where inherent symmetries are present and the case where there
% are no such symmetries and then calculate the SSM coefficients.
%
% If the system is real (the symmetries exist) then only the coefficients up
% to the conjugate center index (in conjugate ordering) are calculated.
%%

Lambda_M = sparse(diag(obj.E.spectrum));
V_M = obj.E.basis;          % Left eigenvectors of the modal subspace
W_M = obj.E.adjointBasis;   % Right eigenvectors of the modal subspace

%% Setup for System with Symmetries

% To make bookkeeping as easy as possible, the eigenvalues and the corresponding
% quantities are sorted such that the $l_i$ complex pairs are in the first $2l_i$
% positions of the linear reduced dynamics, the $l_r$ real eigenvalues are in
% the las $l_r$ positions. This is reverted at the end of computation.
% Sort input according to eigenvalues. Conjugate Evals in the first 2l_i positons, always have to be input in pairs!
% real ones in l_r last pos.
Lambda_M_vector      = diag(Lambda_M);
im_idx= find(imag(Lambda_M_vector) ~= 0); % pos. of imaginary ev pairs in master modal subspace
l_i   = length(im_idx)/2;                 % number of imaginary ev pairs in master modal subspace

r_idx = find(imag(Lambda_M_vector) == 0); % pos. of real eigenvalues in master modal subspace
l_r   = sum(imag(Lambda_M_vector) == 0);  % number of real eigenvalues in master modal subspace

% recast input such that first 2li directions correspond to imag. ev. pairs,
% the last l_r correspond to real eigenvalues.
new_idx    = [im_idx,r_idx];
V_M        = V_M     (:,new_idx);
W_M        = W_M     (:,new_idx);
Lambda_M   = Lambda_M(:,new_idx);
Lambda_M_vector = Lambda_M_vector(new_idx);
%% Conversion of indices and conjugate center index
% To convert between conjugate and reverse lexicographical ordering we construct
% all index sets that do so. If the function computing them is efficient enough
% storing them will no longer be necessary. For now we keep them in the $\texttt{data}$
% field to access them directly.

% conjugate indices
%Z_cci    -  conjugate center indices
%conj2lex -  Set to convert indices from conjugate to rev. lex. ordering
%lex2conj -  Set to convert indices from rev. lex. to conjugate ordering
[multi_input] = conjugate_ordering(order,l_r,l_i);

%% Perparing first order coefficients
% We convert the first order reduced dynamics and SSM-coefficients to conjugate
% ordering to exploit the symmetry in the calculation.
% Sort input in conjugate ordering
idx           = multi_input.revlex2conj{1};
V_M_conj      = V_M(:,idx );    %in conj. ordering
Lambda_M_conj = Lambda_M(:,idx);  %in conj. ordering


V_M_conj(:,1:multi_input.Z_cci(1)); % should be removed
W_0_1(1).coeffs = V_M_conj(:,1:multi_input.Z_cci(1));
R_0_1(1).coeffs = Lambda_M_conj(:,1:multi_input.Z_cci(1));
H{1}  = W_0_1(1).coeffs;




multi_input.H = H;
multi_input.W_M = W_M;
multi_input.Lambda_M_vector = Lambda_M_vector;
multi_input.nl_order = numel(obj.System.F);
multi_input.l = size(Lambda_M,2);

end

function [conjugate] = conjugate_ordering(max_order,l_r,l_i)
% Creates conjugate ordered set indices of multi-indices up to order
% max_order.

% Input:
% max_order   highest order of multi-indices that the conjugate ordering is 
%              wanted for
% l_r number of real eigenvalues in the master subspace
% l_i number of imaignary pairs of eigenvalues in the master subspace

% Output
% for a conjugately ordered set Z, rev. lex ordered set K:
% lex2conj index array for construction.   - K(:,lex2conj) =Z
% conj2lex index array for reconstruction. - Z(:,conj2lex) = K
%     i.e if Z_idx(f) = i, then k_f in K(:,f) is in position i in Z
% cci_k - conjugate center index of the set


Z_cci    = zeros(1,max_order); % conjugate center index array
revlex2conj = cell(1,max_order);  % index sets converting lex set to conj set
conj2revlex = cell(1,max_order);  % index sets converting conj set to lex set

for k = 1:max_order
    I = conjugate_flip(l_i,l_r);
    % Multi_indices in reverse lexicographical ordering
    K = flip(sortrows(nsumk(l_r+2*l_i,k,'nonnegative')).',2);
    Y = K(I,:);
    Exempt      = all(K-Y==0);
    Y(:,Exempt) = [];
    
    %Put m,m_c next to each other
    Z          = [K(:,~Exempt);Y];
    Z          = reshape(Z, size(K,1),[]);
    
    %index out all the combos of m, bar(m) that appear twice
    [~,Z_ia,~] = unique(Z.','rows');
    [~,Z_Ia]   = sort(Z_ia);    
    Z     = Z(:,Z_ia(Z_Ia));
    
    % sort the remaining multi-indices
    idx_1 = 1:2:size(Y,2);
    idx_2 = size(Y,2)-idx_1+1;
    Z     = [Z(:,idx_1),K(:,Exempt),Z(:,idx_2)];
    
    % conjugate center index
    cci_k = sum(Exempt,2) + length(idx_1); 

    % index arrays
    [~,~,conj2lex_k] = intersect(K.',Z.','rows','stable');
    [~,~,lex2conj_k] = intersect(Z.',K.','rows','stable');
    
    conj2revlex{k} = conj2lex_k.';
    revlex2conj{k} = lex2conj_k.';
    Z_cci(k)    = cci_k;
end

conjugate.conj2revlex = conj2revlex;
conjugate.revlex2conj = revlex2conj;
conjugate.Z_cci    = Z_cci;
conjugate.l_i      = l_i;
conjugate.l_r      = l_r;
end

function [idx] = conjugate_flip(l_i,l_r)

% This function computes and index array that flips the conjugate coordinate directions  

idx = reshape(1:2*l_i,2,[]);
idx = reshape(flip(idx),1,[]);
idx = [idx,2*l_i+1:(2*l_i+l_r)];
end