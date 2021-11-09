function AB = tensor_composition(A,B,P,SIZE)
%  TENSOR_COMPOSITION This function computes the composite product of sparse tensor A with the sparse 
% tensors contained in the cell array B in the order specified by the rows of 
% index array P. The returned tensor AB of size SIZE is obtained by summing over 
% the products for each row of indices.
% 
% Let the input $\texttt{A}$  be an $(n+1)-$ dimensional tensor $\mathbf{A}$ 
% and the input $\texttt{B}$ be a $s-$cell array $(s\ge n)$ contain tensors $\mathbf{B}^{(1)},\dots,\mathbf{B}^{(s)}$ 
% and the  input $\texttt{P}$ be an $r\times n$ index matrix $\mathbf{P}=\left[\begin{array}{c}\mathbf{p}_{1}\\\vdots\\\mathbf{p}_{r}\end{array}\right]$. 
% We denote $\mathbf{B}_{\mathbf{p}_j}=\mathbf{B}^{(p_{j1})}\circ\dots\circ\mathbf{B}^{(p_{jn})}$ 
% and compute the tensor composition 
% 
% $$\mathbf{AB} = \sum_{j=1}^{r}\mathbf{A}\cdot\mathbf{B}_{\mathbf{p}_j},$$
% 
% where $\mathbf{A\cdot\mathbf{B} = \mathbf{A}\times_{2}\mathbf{B}^{(1)}\times_3\mathbf{B}^{(2)}\dots\times_{n+1}\mathbf{B}^{(n)}$ 
% is a Tucker-type product. 
AB = sptensor(SIZE);
parfor j_p = 1:size(P,1)
    p = P(j_p,:);
    AB_p = tensor_product(A,B(p));
    AB = AB + AB_p;
end
AB = sptensor(AB);
end

function [AB] = tensor_product(A,B)
%  TENSOR_PRODUCT This function uses the Sandia Tensor Toolbox to compute the product of sparse 
% tensor $\texttt{A}$ with the sparse tensors contained in the cell array $\texttt{B}$.
% 
% % 
% Let the input $\texttt{A}$  be an $(n+1)-$ dimensional tensor $\mathbf{A}$ 
% and the input $\texttt{B}$ be a $n-$cell array contain tensors $\mathbf{B}^{(1)},\dots,\mathbf{B}^{(n)}$ 
% of dimensionality $s_1,\dots,s_n$ respectively. The outer product $\mathbf{B}=\mathbf{B}^{(1)}\circ\dots\circ\mathbf{B}^{(n)}$ 
% has dimensionality $s=s_1+\dots+s_n$. We compute the Tucker-type tensor product 
% $\mathbf{A\cdot\mathbf{B} = \mathbf{A}\times_{2}\mathbf{B}^{(1)}\times_3\mathbf{B}^{(2)}\dots\times_{n+1}\mathbf{B}^{(n)}$ 
% by multiplying the last $n$ modes of $\mathbf{A}$ with the first modes of $\mathbf{B}^{(1)},\dots,\mathbf{B}^{(n)}$ 
% respectively. The resulting tensor has size $n+1-2n+s = s-n+1$. The same is 
% achieved by the following:
AB = A;
for k = 1:length(B)
    AB = sptensor(ttt(AB,B{k},2,1));
end
end