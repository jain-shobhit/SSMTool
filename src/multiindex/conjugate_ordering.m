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