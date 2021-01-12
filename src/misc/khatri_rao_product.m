function AB = khatri_rao_product(A,B)
%% 
% This function computes the column-wise khatri-rao product of matrices A, B

[m_A, n_A] = size(A);
[m_B, n_B] = size(B);
if n_A ~= n_B
    error('Both matrices must have the same number of columns')
end
if issparse(A) || issparse(B)
    AB = sparse(m_A*m_B,n_A);
else
    AB = zeros(m_A*m_B,n_A);
end
for j = 1:n_A
    AB(:,j) = kron(A(:,j),B(:,j));
end
end