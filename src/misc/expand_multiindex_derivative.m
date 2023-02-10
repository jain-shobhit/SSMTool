function S = expand_multiindex_derivative(M,P)
%  EXPAND_MULTIINDEX This function returns the derivative of multi-index expansion $\mathbf{S}$ for the multi-index
% pair $(\mathbf{C},\mathbf{I})$ over the column vectors in the matrix $\mathbf{P}=[\mathbf{p}_1,\dots,\mathbf{p}_{n_t}]\in
% \mathbb{R}^{n\times n_t}$, where $\mathbf{C},\mathbf{I}$ are fields of
% the data structure M.
%
% $\mathbf{s}_j = D(\mathbf{C} (\mathbf{p}_j)^{\mathbf{I}})/Dp,\quad j=1,\dots,n_t$,
%
% $$\mathbf{S}=[\mathbf{s}_1,\dots,\mathbf{s}_{n_t}]$$ (in 3D array)
if isempty(M.ind)
    S = 0;
else
    nt = size(P,2);
    C = M.coeffs;
    I = M.ind;
    n_I = size(I,1);
    dimSys = size(C,1);
    dimSSM = size(I,2);
    S = zeros(dimSys,dimSSM,nt);
    for i=1:n_I
        coeff = C(:,i);
        ind = I(i,:)';
        % find nonzero exponents
        expind = find(ind);
        pind   = P(expind,:);
        pind(pind==0) = eps; % handle divided by zero
        s = prod(pind.^ind(expind),1);    
        s = (ind(expind)./pind).*s; % #p X nt
        s = s.'; s = s(:); s = s.';
        s = coeff.*s; s = reshape(full(s),[dimSys,nt,numel(expind)]); s = permute(s,[1,3,2]);
        S(:,expind,:) = S(:,expind,:)+ s;
    end
end
end

