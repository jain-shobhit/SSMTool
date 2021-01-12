function S = expand_multiindex(M,P)
%  EXPAND_MULTIINDEX This function returns the multi-index expansion $\mathbf{S}$ for the multi-index
% pair $(\mathbf{C},\mathbf{I})$ over the column vectors in the matrix $\mathbf{P}=[\mathbf{p}_1,\dots,\mathbf{p}_{n_t}]\in
% \mathbb{R}^{n\times n_t}$, where $\mathbf{C},\mathbf{I}$ are fields of
% the data structure M.
%
% $\mathbf{s}_j = \mathbf{C} (\mathbf{p}_j)^{\mathbf{I}},\quad j=1,\dots,n_t$,
%
% $$\mathbf{S}=[\mathbf{s}_1,\dots,\mathbf{s}_{n_t}]$$
if isempty(M.ind)
    S = 0;
else
    nt = size(P,2);
    C = M.coeffs;
    I = M.ind;
    n_I = size(I,1);
    % evaluate monomials
    P_I = prod(kron(P.',ones(n_I,1)) .^ repmat(I,nt,1),2);
    % evaluate polynomials
    S = C*reshape(P_I,[],nt);
end
end

% % loop version
% function S = expand_multiindex(M,P)
% %  EXPAND_MULTIINDEX This function returns the multi-index expansion $\mathbf{S}$ for the multi-index
% % pair $(\mathbf{C},\mathbf{I})$ over the column vectors in the matrix $\mathbf{P}=[\mathbf{p}_1,\dots,\mathbf{p}_{n_t}]\in
% % \mathbb{R}^{n\times n_t}$, where $\mathbf{C},\mathbf{I}$ are fields of
% % the data structure M.
% %
% % $\mathbf{s}_j = \mathbf{C} (\mathbf{p}_j)^{\mathbf{I}},\quad j=1,\dots,n_t$,
% %
% % $$\mathbf{S}=[\mathbf{s}_1,\dots,\mathbf{s}_{n_t}]$$
% if isempty(M.ind)
%     S = 0;
% else
%     n_I = size(I,1);
%     nt = size(P,2);
%     C = M.coeffs;
%     I = M.ind;
%     
%     Pprod=zeros(n_I,nt);
%     for j = 1:nt
%         Pj = repmat(P(:,j).',[n_I,1]);
%         Pprod(:,j) = prod(Pj.^I,2);
%     end
%     S = C*Pprod;
% end
% end