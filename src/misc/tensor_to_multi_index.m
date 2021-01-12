function M = tensor_to_multi_index(T)
% Convert polynomials in tensor format to multiindex format. Here T is a
% single tensor or a cell array of tensors. For each tensor, the first mode
% gives the index of polynomial vector, and the rest modes correspond to
% factors of polymonials, e.g. x^3=x*x*x.
% The result is return by two matrices, coeffs and ind, which represents
% coefficient and multiindex respesctively. Specifically, each column in
% coeffs gives the coefficient vector for a polynomial whose indices are
% stored as a row vector in ind matrix.

M.coeffs = [];
M.ind = [];
if ~iscell(T)
    [M.coeffs,M.ind] = get_multiindices(T);
else
    
for j = 1:length(T)    
    [C_j,I_j] = get_multiindices(T{j});
    M.coeffs = cat(2,M.coeffs,C_j);
    M.ind = cat(1,M.ind,I_j);    
end

end
end
%%
function [C,I] = get_multiindices(T)
%% 
% This function returns the multi-index coefficients and the corresponding exponents 
% for the given set of multi-dimensional arrays 
% 
% Outputs: 
%% 
% * $\mathbf{C} \in \mathbb{R}^{ n\times n_I}$: Array containing multi-index 
% coefficients
% * $\mathbf{I} \in \mathbb{R}^{ n_I\times r}$:  Array containing multi-index 
% exponents
%% 
% Inputs:
%% 
% * $\mathbf{T} \in \mathbb{R}^{ n\times r\times \dots\times r}$: $(d+1)$-dimensional 
% array in the form of sptensor object (cf. tensor toolbox)
%% 
% Other variables:
%% 
% * $d\in\mathbb{N}$: degree of tensor polynomial (only uniform single 
% degree supported currently)
% * $r\in\mathbb{N}$: The numer of variables in the multivariate polynomial     
%% 
% Assuming the tensor array $\mathbf{T}$ acts on the vector $\mathbf{q}\in\mathbb{R}^{r}$ 
% along the dimensions $2,\dots,d+1$. We have the following relationship between 
% the input and output
% 
% $$\sum_{j_1,\dots,j_d=1}^r T_{ij_1\dots j_d} q_{j_1} q_{j_2} \dots q_{j_d} 
% = \sum_{\mathbf{p}\in\mathbb{N}_0^{r},~|\mathbf{p}|=d}C_{i,\mathbf{p}}\mathbf{q}^{\mathbf{p}},\quad 
% \mathbf{q}\in\mathbb{R}^{r},\quad i\in \{1,\dots,n\},\quad d = \textrm{dim}(\texttt{T}) 
% - 1$$ 
%
% where $\mathbf{q}^{\mathbf{p}}=q_1^{p_1}\cdots q_r^{p_r}$
% 
% or
%% $\mathbf{T}\mathbf{q}^{\otimes d} = \mathbf{C}\mathbf{q}^{\mathbf{I}}$,
% 
% where the different arrays are organized in the following form:
% 
% $\mathbf{C} =  \left[\begin{array}{c}\mathbf{c}_1\\\vdots\\\mathbf{c}_{n}\end{array}\right]$, 
% where the $i^{\mathrm{th}}$ row represents the coefficients $\mathbf{c}_i = 
% [C_{i,\mathbf{p}_1},\dots,C_{i,\mathbf{p}_{n_I}}]$,
% 
% $\mathbf{I} = \left[\begin{array}{c}\mathbf{p}_1\\\vdots\\\mathbf{p}_{n_I}\end{array}\right]$ 
% is a matrix whose rows contain the multi-indices $\mathbf{p}\in\mathbb{N}_0^{r}$ 
% respresented by the tensor $\mathbf{T}$, 
% 
% $\mathbf{q}^{\mathbf{I}} = \left[\begin{array}{c}\mathbf{q}^{\mathbf{p}_1}\\\vdots\\\mathbf{q}^{\mathbf{p}_{n_I}}\end{array}\right]$.
n = size(T,1);
r = size(T,2);
[SUBS, VALS] = sparsify(T);
if ~isempty(SUBS)
    %% %%% Old version %%%
    %{
    multi_indices = sub2multiind(SUBS(:,2:end), r);
    [I,~,IC] = unique(multi_indices,'rows');
    
    n_I = size(I,1);
    C = sparse(SUBS(:,1),IC,VALS,n,n_I);
    %}
    
    %% %%% New version %%%
    % Instead of creating the multi-indices for all subscripts and then
    % finding unique multi-indices this version finds all unique subscripts
    % up to permutations 
    % eg. [1,0,2] and [1,2,0] correspond to the same multi-index and thus
    % the multi-indices are created for only one of them.
    
    % The main advantage of this method comes from the fact that at order k
    % there are only k subscripts for a multi-index, but the corresponding 
    % multi-indices have the dimension of the space the tensor T maps into 
    % which may be very high. Performing operations such as unique on multi-indices is thus
    % much more expensive than finding the subscripts that those unique
    % multi-indices correspond to.
    % {
    
    % sort each subscript such that subscripts that are the same up to
    % permutations are identical.
    SUBS_un = sort(SUBS(:,2:end),2);
    % sort subscripts such that identical subscripts are placed in adjacent rows    
    [SUBS_un, sort_idx] = sortrows(SUBS_un);
    
    % by finding the subscripts that are different from their neighbour we
    % find the position of all unique subscripts. Note that on top of the
    % array SUBS_un we insert a row of zeros to get the right indices since
    % diff(r) outputs an array with one row less than r
    subs_un_idx = any(diff([zeros(1,size(SUBS,2)-1); SUBS_un]),2);

    % The array multi_indices only contains
    % unique multi-indices the resulting multi-indices are in reverse order compared to the
    % output of the old version.
    multi_indices = sub2multiind(SUBS_un(subs_un_idx,:), r);
    % same ordering of multi-indices as in the old version
    I = flip(multi_indices,1); 
    n_I = size(I,1);
    
    % IC contains the position of the multi-index in multi_indices
    % that a subscript in SUBS_un corresponds to for all subscripts
    IC = cumsum(subs_un_idx);
    
    % Since the ordering of the subscripts in IC is not the same as in SUBS
    % we have to construct an index array that reverts the sorting that has
    % been previously applied
    sort_idx_rev     = zeros(1,length(sort_idx));
    sort_idx_rev(sort_idx) = 1:length(sort_idx);
    
    % Using the reverse index array we put the indices in IC in the
    % position of the subscript they correspond to in SUBS
    IC = IC(sort_idx_rev);
    
    % Coefficients are read out in the multi-index format 
    C = sparse(SUBS(:,1),IC,VALS,n,n_I);
    % Also put coefficients in same ordering as the output of the old version
    C = flip(C,2); 
    %}
else 
    C = [];
    I = [];
end
end
%% 
% 
function [subs,vals] = sparsify(T)
if isa(T,'sptensor')
    subs = T.subs;
    vals = T.vals;
else
    SIZE = size(T);
    subs = tt_ind2sub(SIZE, (1:prod(SIZE))');
    vals = T(:);
end

end