function[Hk] = composition_coefficients_of_powerseries(W_0,H,multi_input)
%  COMPOSITIONCOEFFICIENTS This function computes the coefficients for the composition of a multivariate
%
% polynomial power series.
%
% Take the expansion of row $i$ of a vector valued function $(\mathbf{S(p)}
% )_i= \sum_{\mathbf{u}}S_{\mathbf{u},i}\mathbf{p^u}$ .
%
% If this power series is multiplied with itself $s$ times, it can be written
% like $((\mathbf{S(p)} )_i)^s= \bigg( \ \sum_{\mathbf{u}}S_{\mathbf{u},i}\mathbf{p^u}\bigg)^s=
% \sum_{\mathbf{h}}H_{i,s,\mathbf{h}}\mathbf{p^h}$
%
% This function computes $H_{i,s,\mathbf{k}_f} = \frac{s}{k_{f,j}} \sum_{ \mathbf{u
% < k}_f} u_j S_{i,\mathbf{u}} H_{i,s-1,\mathbf{k}_f-\mathbf{u}} $ for all multi-index
% vectors $\mathbf{k}_f$ of a given order $k$.
Z_cci = multi_input.Z_cci;
N = multi_input.N;
K = multi_input.K;
revlex2conj = multi_input.revlex2conj;
l = size(K,1);
k = sum(K(:,1),1);

string = 'conjugate';
z_k    = Z_cci(k);

Hk = zeros(N,z_k,k);
[kj,Ikj] = min(K +(k+1)*(K ==0)); %find nonzero minima of all multi-indices of order k
%%
% Loop over all orders of multi-indices that are potentially nonzero and do
% not require the knowledge of the order $k$ coefficients of the SSM parametrisation.
for ord = 1:k-1
    %%
    % Create all multi-indices of order $k-u$
    g       = flip(sortrows(nsumk(l,k-ord,'nonnegative')).',2);
    
    g       = g(:,revlex2conj{k-ord});
    
    [gmink,I_k,I_g] = gminusk(g,K);
    %%
    % Now all the positions of the order $m$ multi-indices are calculated. Assume
    % the vector $\texttt{gmink}$ has $u_\alpha$ columns.
    I_gmk     = multi_index_2_ordering(gmink,string,revlex2conj);
    %%
    % We read out all the values $u_j$ positions of minimal entries
    pos     = Ikj(I_k);                 %find row for every vector in mv
    sz      = size(I_k,2);
    poslin  = sub2ind([l,sz],pos,1:sz); %create linear index array
    % read out minimal entries
    gmk_j     = gmink(poslin);            %has to be a row vector
    %%
    % One thing that is redundant is all the calculations where $u_j$ is zero. The
    % corresponding values could be ignored, by taking the indices corresponding to
    % $u_j$= 0 out of calculations .
    %
    % Consequently, the coefficients $u_jS_{i,\mathbf{u}} H_{i,s-1,\mathbf{k}_f-\mathbf{u}}
    % $  are calculated. In the array $\texttt{Hk}$ the first row corresponds to the
    % phase space dimension, the third one to the scalar exponents.
    %
    % $\texttt{H}_{k-u}(:,\texttt{I}_{g},:)$ is a matrix in $\mathbb{C}^{2n \times
    % \texttt{sz} \times k-1}$. Contribution of coefficients where $s=1$ explicitly
    % depend on the SSM coefficients of the order that is calculated and will be added
    % in the end after the highest order coefficients of the SSM parametrisation have
    % been computed,but fortunately are not needed for the calculations of themselves.
    % Finally $\texttt{S\{ord\}}(:,\texttt{gmk\_I}) \in \mathbb{C}^{2n \times \texttt{sz}
    % }$.
    
    %%
    % Now, since we want to make use of the symmetry in the SSM coefficients, we
    % split all contributions such that in conjugate multi-index ordering, all of
    % the indices at order $\texttt{ord}$ that are smaller or equal than $\texttt{Z(ord)}$
    % are read out normally, and the ones that are bigger are changed to their conjugate
    % counterpart.
    coeff       = zeros(N,sz,k-ord);
    
    z_o   = nchoosek(ord+l-1,l-1);
    z_ko  = nchoosek(k-ord+l-1,l-1);
    z_cci_o = Z_cci(ord);
    z_cci_ko = Z_cci(k-ord);
    
    %determine which coefficients correspond to multi-indices bigger
    %than the conjugate center index at order o and k-o.
    
    I_g_a = I_g <= z_cci_ko;
    I_gmk_a = I_gmk <= z_cci_o;
    % where g_I <= z_ko_2, and gmk_I <= z_k_2
    aa = I_g_a & I_gmk_a;
    % where g_I > z_ko_2, but gmk_I <= z_k_2
    ab = ~I_g_a & I_gmk_a;
    % where g_I <= z_ko_2, but gmk_I > z_k_2
    ba = I_g_a & ~I_gmk_a;
    % where g_I > z_ko_2, and gmk_I > z_k_2
    bb = ~ I_g_a & ~I_gmk_a;
    
    %read out contributions using symmetries
    coeff(:,aa,:)   = W_0{ord}(:,I_gmk(aa)).* gmk_j(aa)             .* H{k-ord}(:,I_g(aa),:);
    coeff(:,ab,:)   = W_0{ord}(:,I_gmk(ab)).* gmk_j(ab)             .* conj(H{k-ord}(:,z_ko-I_g(ab)+1,:));
    coeff(:,ba,:)   = conj(W_0{ord}(:,z_o-I_gmk(ba)+1)).* gmk_j(ba) .* H{k-ord}(:,I_g(ba),:);
    coeff(:,bb,:)   = conj(W_0{ord}(:,z_o-I_gmk(bb)+1)).* gmk_j(bb) .* conj(H{k-ord}(:,z_ko-I_g(bb)+1,:));
    
    %%
    % $\texttt{coeff} \in \mathbb{C}^{2n \times \texttt{sz} \times k-1}$. What is
    % left now is to add up all the parts corresponding to a specific multi-index
    % of order  $k$ and multiply with the $s$ the slice corresponds to. We want to
    % add all the columns corresponding to where $\texttt{I}_{k}$ has the same values
    % (the columns contain coefficients that contribute to the same order $k$ multi-index).
    for f = unique(I_k)
        Hk(:,f,2:(k-ord+1)) = Hk(:,f,2:(k-ord+1)) + sum(coeff(:,I_k==f,:),2);
    end
end
%multiplication with s and kj
s  = 1:k;
%match dimensions for pointiwse multiplication
Hk = permute(Hk,[2,3,1]) .* ((1./kj).' * s);
Hk = permute(Hk,[3,1,2]);
end