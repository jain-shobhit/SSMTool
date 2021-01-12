function [T] = multi_index_to_tensor(C,I)
%  MULTI_INDEX_TO_TENSOR This function returns the multi-dimensional array for given multi-index coefficients 
% and the corresponding exponents
% 
% Inputs: 
%% 
% * $\mathbf{C} \in \mathbb{R}^{ n\times n_I}$: Array containing multi-index 
% coefficients
% * $\mathbf{I} \in \mathbb{R}^{ n_I\times r}$:  Array containing multi-index 
% exponents
%% 
% Outputs:
%% 
% * $\mathbf{T} \in \mathbb{R}^{ n\times r\times \dots\times r}$: $(d+1)$-dimensional 
% arraty in the form of sptensor object (cf. tensor toolbox)
%% 
% Other variables:
%% 
% * $d\in\mathbb{N}:$ degree of tensor polynomial (currently only uniform single 
% degree supported currently
% * $r\in\mathbb{N}:$ The numer of variables in the multivariate polynomial     
%% 
% Assuming the tensor array $\mathbf{T}$ acts on the vector $\mathbf{q}\in\mathbb{R}^{r}$ 
% along the dimensions $2,\dots,d+1$. We have the following relationship between 
% the input and output
% 
% $$\sum_{j_1,\dots,j_d=1}^r T_{ij_1\dots j_d} q_{j_1} q_{j_2} \dots q_{j_d} 
% = \sum_{\mathbf{p}\in\mathbb{N}_0^{r},~|\mathbf{p}|=d}C_{i,\mathbf{p}}\mathbf{q}^{\mathbf{p}},\quad 
% \mathbf{q}\in\mathbb{R}^{r},\quadi\in \{1,\dots,n\},\quad d = \textrm{dim}(\texttt{T}) 
% - 1$$ 
% 
% or
%% $\mathbf{T}\mathbf{q}^{\otimes d} = \mathbf{C}\mathbf{q}^{\mathbf{I}}$,
% % 
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
n = size(C,1);
d = sum(I(1,:));
r = size(I,2);
nSubs = size(I,1);
SIZE = [n r*ones(1, d)];
SUBS = zeros(nSubs,d);
for j = 1:nSubs
    p = I(j,:);
    SUBS(j,:) = multi_index_to_sub(p);
end
[i,j,VALS] = find(sparse(C));
T = sptensor([i, SUBS(j,:)],VALS,SIZE);
end
%%
function subs = multi_index_to_sub(p)
    d = sum(p);
    subs = zeros(1,d);
    [~, i, s] = find(p);
    m = 1;
    for k = 1:d
        subs(k) = i(m);
        if s(m) == 1
            m = m + 1;
        else
            s(m) = s(m) - 1;
        end
    end
end