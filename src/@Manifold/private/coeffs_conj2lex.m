function [W_0_lex,R_0_lex] = coeffs_conj2lex(multi_input,order,W_0,R_0)
%%
% This function reconstructs the full SSM-coefficients and reduced dyn.
% in lexicographical ordering from the conjugate ordering coefficients
%%
l   = multi_input.l; % SSM dimension and partition into real and imaginary subspaces
l_i = multi_input.l_i;
l_r = multi_input.l_r;

Z     = number_of_multis(l,order); % Total number of size l multi-indices for every order
Z_cci = multi_input.Z_cci; % conjugate center indices, denoted as 

col_idx = num2cell(Z-Z_cci);
row_idx = mat2cell(repmat(conjugate_flip(l_i,l_r),1,order),1,l*ones(1,order));

% Index array to convert conjugate to lexicographical ordering
conj2revlex = multi_input.conj2revlex;
conj2lex    = cellfun(@flip , conj2revlex, num2cell(2*ones(1,order)), 'UniformOutput',false);  %to convert to lexicographic ordering

% convert to cell to apply cellfun
W_0 = num2cell(W_0);
R_0 = num2cell(R_0);

[type{1:order}] = deal('TaylorCoeffs');
W_0_lex     = cellfun(@coeffs_conj2full, W_0,cell(1,order), col_idx,conj2lex,type, 'UniformOutput',false);
R_0_lex     = cellfun(@coeffs_conj2full, R_0,row_idx , col_idx,conj2lex, type, 'UniformOutput',false);

%convert back to struct array

W_0_lex = [W_0_lex{:}];
R_0_lex = [R_0_lex{:}];


%[type{1:order}] = deal('CompCoeffs');
%H	        = cellfun(@coeffs_conj2full, multi_input.H ,cell(1,order), col_idx,conj2lex, type, 'UniformOutput',false);
%save('H.mat','H');
end

function [Z] = number_of_multis(l,max_order)

%Returns the number of size l multi-indices for all orders up to order max_order

Z = zeros(1,max_order);

for j = 1:max_order
    Z(j) = nchoosek(j+l-1,l-1);
end
end

function [idx] = conjugate_flip(l_i,l_r)

% This function computes and index array that flips the conjugate coordinate directions  

idx = reshape(1:2*l_i,2,[]);
idx = reshape(flip(idx),1,[]);
idx = [idx,2*l_i+1:(2*l_i+l_r)];
end