function [L_k] = mixed_terms_in_SR(W_0,R_0,multi_input)
%% MIXEDTERMS Mixed Terms in the Invariance Equation
% This function computes all the mixed terms of the invariance equation at order
% $k$ for a $l$ dimensional SSM.
%
% The mixed terms on the right hand side of the invariance equation are given
% by all contributions of SSM-coefficients that do not have linear order or order
% $k$.
Z_cci       = multi_input.Z_cci;
N           = multi_input.N;
K           = multi_input.K;

revlex2conj = multi_input.revlex2conj;
l           = size(K,1);
k           = sum(K(:,1),1);
%%
% If the system is real we need all lower order conjugate center indices which
% are stored in $\texttt{Z}$. If the system is not real, only the number of multi-indices
% at order $k$ is needed and thus $\texttt{Z}$ is empty.
string = 'conjugate';
z_k    = Z_cci(k);

L_k = sparse(N,z_k);
for u = 2:k-1
    % Anti-lex. ordered multi-indices of order u
    %%
    % For every one of these order $u$ multi-indices, the coefficients $u_j S_{i,\mathbf{u}}R_{j,\mathbf{g}}$
    % are calculated.
    multis = flip(sortrows(nsumk(l,u,'nonnegative')).',2);
    
    multis = multis(:,revlex2conj{u});
    z_u_2 = Z_cci(u);
    z_u = nchoosek(u+l-1,l-1);
    z_g =nchoosek(k-u+1+l-1,l-1);
    
    
    for j  = 1:l
        ej = (1:l == j).';
        
        %find j for using the symmetry in reduced dynamics
        if mod(j,2) == 0
            jbar = j-1;
        else
            jbar = j+1;
        end
        %%
        % This loop goes over all of the multi-index vectors at order $u$, and in every
        % loop all coefficients of the mixed terms corresponding to $\mathbf{g} = \mathbf{k}_f
        % + \mathbf{e}_j - \mathbf{u}_{index}$are determined for all $f = 1,...z_k$.
        for index = 1:size(multis,2)
            if  multis(j,index) ~= 0
                [g,k_idx,~] = gminusk(multis(:,index),K+ej);
                if ~isempty (g)
                    g_idx        =multi_index_2_ordering(g,string,revlex2conj);
                    
                    % Make use of the Symmetry
                    
                    g_I_a = g_idx <= Z_cci(k+1-u);
                    term  = sparse(N,size(k_idx,2));
                    
                    if index <= z_u_2
                        term(:,g_I_a)  = multis(j,index) * W_0{u}(:,index) .* R_0{(k-u+1)}(j,g_idx(g_I_a));
                        term(:,~g_I_a) = multis(j,index) * W_0{u}(:,index) .* conj(R_0{(k-u+1)}(jbar,z_g-g_idx(~g_I_a)+1));
                    else
                        
                        term(:,g_I_a)  = multis(j,index) * conj(W_0{u}(:,z_u-index+1)) .* R_0{(k-u+1)}(j,g_idx(g_I_a));
                        term(:,~g_I_a) = multis(j,index) * conj(W_0{u}(:,z_u-index+1)) .* conj(R_0{(k-u+1)}(jbar,z_g-g_idx(~g_I_a)+1));
                    end
                    
                    % Sum contributions
                    for f = unique(k_idx)
                        L_k(:,f) = L_k(:,f) - sum(term(:,k_idx==f),2);
                    end
                end
            end
            
        end
    end
end