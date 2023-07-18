function [res] = coeffs_mixed_terms(k,m,W,R,field,mix)
%%
% This function computes all the mixed terms of the invariance equation at order
% $k$ for a $l$ dimensional SSM. It is used for nonautonomous computation
% in revlex ordering and for the autonomous SSM computation in conjugate
% ordering.

% In the case where the DS is complex then the autonomous computation also
% uses the revlex ordering.
%%

N   = field.N;
l   = field.l;


switch field.ordering
    case 'revlex'
        z_k = nchoosek(k+l-1,l-1);
        res = zeros(N,z_k);
        switch mix
            % strictly doesnt need distinction since u is the same but due
            % to different reasons
            case 'R1'
                
                u = k + 1 - m + 1; %index into non-autonomous, starts at 0
                
            case 'W1'
                
                u = k + 1 -(m - 1); %index into autonomous, starts at 1
            case 'aut'
                
                u = k - m + 1; % In this case there is no zeroth order term
        end
        
        if nnz(W(m).coeffs)>0 && nnz(R(u).coeffs)>0
            
            if l >1
                K   = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
                if isfield(W,'ind') % nonautonomous computation
                    M = W(m).ind; %order m or m+1 (W1)multi-indices
                    U = R(u).ind; % order u or u+1 (R1) multi- indices
                else % autonomous computation
                    M = flip(sortrows(nsumk(l,m,'nonnegative')).',2); %order m multi-indices
                    U = flip(sortrows(nsumk(l,u,'nonnegative')).',2); %order u multi-indices
                end
            else
                K = k;
                if isfield(W,'ind') % nonautonomous computation
                    M = W(m).ind; %order m or m+1 (W1)multi-indices
                    U = R(u).ind; % order u or u+1 (R1) multi- indices
                else % autonomous computation
                    M = m; %order m multi-indices
                    U = u; %order u multi-indices
                end
            end

            One = eye(l);
            
            [MU,I_m,I_u] = multi_addition(M,U);
            
            for j = 1:l
                
                %copy for this run
                MU_One = MU - One(:,j);
                I_m_One = I_m;
                I_u_One = I_u;
                
                %remove negative elements
                [~, neg] = find (MU_One<0);
                MU_One(:,neg)= [];
                I_m_One(:,neg)= [];
                I_u_One(:,neg)= [];
                
                %Assemble coefficients
                Wm = W(m).coeffs(:,I_m_One);
                mRu = R(u).coeffs(j,I_u_One).* M(j,I_m_One);
                WR  = Wm .* mRu;
                
                %find position of multi-indices in the set K
                [~,loc_k] = ismember(MU_One.',K.','rows');
                
                for f = 1:z_k
                    if any(loc_k==f)
                        res(:,f) = res(:,f) + sum( WR(:,loc_k==f),2);
                    end
                end
            end
            
        end
        
        
    case 'conjugate'
        
        % The mixed terms on the right hand side of the invariance equation are given
        % by all contributions of SSM-coefficients that do not have linear order or order
        % $k$.
        Z_cci       = field.Z_cci;
        K           = field.K;
        revlex2conj = field.revlex2conj;
        k           = sum(K(:,1));
        %%
        % If the system is real we need all lower order conjugate center indices which
        % are stored in $\texttt{Z}$. If the system is not real, only the number of multi-indices
        % at order $k$ is needed and thus $\texttt{Z}$ is empty.
        string = 'conjugate';
        z_k    = Z_cci(k);
        
        res = sparse(N,z_k);
        % Anti-lex. ordered multi-indices of order u
        %%
        % For every one of these order $u$ multi-indices, the coefficients $u_j S_{i,\mathbf{u}}R_{j,\mathbf{g}}$
        % are calculated.
        M = flip(sortrows(nsumk(l,m,'nonnegative')).',2);
        
        M = M(:,revlex2conj{m});
        z_u_2 = Z_cci(m);
        z_u = nchoosek(m+l-1,l-1);
        z_g = nchoosek(k-m+1+l-1,l-1);
        
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
            for index = 1:size(M,2)
                if  M(j,index) ~= 0
                    [g,k_idx,~] = multi_subtraction(K+ej,M(:,index),'Parametrised');
                    if ~isempty (g)
                        g_idx        =multi_index_2_ordering(g,string,revlex2conj);
                        
                        % Make use of the Symmetry
                        
                        g_I_a = g_idx <= Z_cci(k+1-m);
                        term  = sparse(N,size(k_idx,2));
                        
                        if index <= z_u_2
                            term(:,g_I_a)  = M(j,index) * W(m).coeffs(:,index) .* R((k-m+1)).coeffs(j,g_idx(g_I_a));
                            term(:,~g_I_a) = M(j,index) * W(m).coeffs(:,index) .* conj(R((k-m+1)).coeffs(jbar,z_g-g_idx(~g_I_a)+1));
                        else
                            
                            term(:,g_I_a)  = M(j,index) * conj(W(m).coeffs(:,z_u-index+1)) .* R((k-m+1)).coeffs(j,g_idx(g_I_a));
                            term(:,~g_I_a) = M(j,index) * conj(W(m).coeffs(:,z_u-index+1)) .* conj(R((k-m+1)).coeffs(jbar,z_g-g_idx(~g_I_a)+1));
                        end
                        % Sum contributions
                        for f = unique(k_idx)
                            res(:,f) = res(:,f) + sum(term(:,k_idx==f),2);
                        end
                    end
                    
                end
            end
        end
end