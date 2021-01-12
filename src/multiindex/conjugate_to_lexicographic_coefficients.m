function [W_0_lex,R_0_lex] = conjugate_to_lexicographic_coefficients(multi_input,order,W_0,R_0)

%%
%If the system is real, the following code reconstructs the full SSM-coefficients
% in lexicographical ordering.

l   = multi_input.l; % SSM dimension and partition into real and imaginary subspaces
l_i = multi_input.l_i;
l_r = multi_input.l_r;

Z     = number_of_multis(l,order); % Total number of size l multi-indices for every order
Z_cci = multi_input.Z_cci; % conjugate center indices, denoted as 

idx_1 = num2cell(Z-Z_cci);
idx_2 = mat2cell(repmat(conjugate_flip(l_i,l_r),1,order),1,l*ones(1,order));

% Index array to convert conjugate to lexicographical ordering
conj2revlex = multi_input.conj2revlex;
conj2lex    = cellfun(@flip , conj2revlex, num2cell(2*ones(1,order)), 'UniformOutput',false);  %to convert to lexicographic ordering

% All coefficients in lexicographic ordering of the corresponding
% multi-indices
W_0_lex     = cellfun(@restore_full_coeff, W_0,cell(1,order), idx_1,conj2lex, 'UniformOutput',false);
R_0_lex     = cellfun(@restore_full_coeff, R_0,idx_2 , idx_1,conj2lex, 'UniformOutput',false);

end