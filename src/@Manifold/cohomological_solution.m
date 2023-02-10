function [W_0i, R_0i,multi_input] = cohomological_solution(obj, i,  W_0, R_0,multi_input)
%% COHOMOLOGICAL_SOLUTION Solution of cohomological equations at order i
switch obj.Options.notation
    case 'tensor'
        % The function computes the SSM where we solve the invariance equation
        %
        % $$\mathbf{B}D\mathbf{S}\mathbf{R}=\mathbf{A}\mathbf{S}+\mathbf{F}\circ\mathbf{S}$$
        %
        % of the dynamical system
        %
        % $\mathbf{B}\dot{\mathbf{z}} = \mathbf{A}\mathbf{z} + \mathbf{F(\mathbf{z})}$.
        %
        % The SSM is expressed in terms of the expansion
        %
        % $$\mathbf{S}(\mathbf{p})=\sum_{i=1}^{\Gamma_{S}}\mathbf{S}_{i}\mathbf{p}^{\otimes
        % i},$$
        %
        % where $\mathbf{p}\in\mathbb{C}^m$ are parameterization coordinates of the
        % $m-$dimensional SSM.
        %
        % The coefficients at different orders are collected in a cell array $\texttt{S}$
        % , where $\texttt{S\{i\}}$ gives the coefficients at order $\texttt{i}$, i.e.
        % $\mathbf{S}_{i}$. These are obtained by solving for $\mathbf{S}_{i}$ in the
        % following equation
        %
        % $$\mathbf{B}\mathbf{S}_{i}\mathbf{\Lambda}_{\mathcal{M},i}-\mathbf{A}\mathbf{S}_{i}=\underbrace{\sum_{j=2}^{i}\mathbf{F}_{j}\sum_{|\mathbf{p}|=i}\mathbf{S}_{p_{1}}\otimes\dots\otimes\mathbf{S}_{p_{j}}-\mathbf{B}\sum_{j=2}^{i-1}\mathbf{S}_{j}\sum_{|\mathbf{p}|=1}\mathbf{R}_{i+1-j}^{p_{1}}\otimes\dots\otimes\mathbf{R}_{i+1-j}^{p_{j}}}_{\mathbf{L}_{i}}-\mathbf{B}\mathbf{S}_{1}\mathbf{R}_{i}$$
        %
        % where
        %
        % $$\mathbf{\Lambda}_{\mathcal{M},i}:=\sum_{|\mathbf{p}|=1}\mathbf{\Lambda}_{\mathcal{M}}^{p_{1}}\otimes\dots\otimes\mathbf{\Lambda}_{\mathcal{M}}^{p_{i}}$$
        %
        % and $\mathbf{\Lambda}_{\mathcal{M}}$ is a diagonal ($m\times m$) matrix containing
        % the eigenvalues of the master modal subspace $\mathcal{M}$. The above equation
        % in the vectorized notation is given by
        %
        % $$\underbrace{\left\left[\left(\mathbf{\Lambda}_{\mathcal{M},i}^{\top}\otimes\mathbf{B}\right)-\left(\mathbf{I}_{m^{i}}\otimes\mathbf{A}\right)\right]}_{{\mathcal{C}}_{i}}\mathfrak{vec}\left(\mathbf{S}_{i}\right)=\mathfrak{vec}\left(\mathbf{L}_{i}\right)-\left(\mathbf{I}_{m^{i}}\otimes\mathbf{B}\mathbf{S}_{1}\right)\mathfrak{vec}\left(\mathbf{R}_{i}\right)$$
        %
        % Here $\mathfrak{vec}\left(\mathbf{S}_{i}\right)$ just stands for the vectorization
        % operator in MATLAB obtained by the command $\texttt{S\{i\}(:)}$.
        %
        Lambda_M = obj.E.spectrum;
        A = obj.System.A; % A matrix
        B = obj.System.B; % B matrix
        W_M = obj.E.adjointBasis; % Right eigenvectors of the modal subspace
        V_M = obj.E.basis; % Left eigenvectors of the modal subspace
        N = obj.dimSystem; % Full system dimensionality in first-order form
        F = obj.System.F; % Full system Nonlinearity coefficients at different orders
        m = length(Lambda_M); % dim(M): M is the master modal subspace
        N_i = N*m^i; % number of unknown SSM coefficients in the tensor notation at order i
        ref = min(abs(Lambda_M));
        if ref<1e-10; ref = max(abs(Lambda_M)); end
        abstol = obj.Options.reltol * ref;
                
        %% Assemble the coefficient matrix of SSM
        % Obtaining $\mathbf{\Lambda}_{\mathcal{M},i}:=\sum_{|\mathbf{p}|=1}\mathbf{\Lambda}_{\mathcal{M}}^{p_{1}}\otimes\dots\otimes\mathbf{\Lambda}_{\mathcal{M}}^{p_{i}}$
        %
        % We assemble it as
        %
        % $$\mathbf{\Lambda}_{\mathcal{M},i}=\sum_{j=1}^{m}\mathbf{I}_{m}\otimes\dots\otimes\mathbf{I}_m\otimes\mathbf{\Lambda}_{\mathcal{M}}\otimes\mathbf{I}_m\otimes\dots\otimes\mathbf{I}_{m}\,,$$
        %
        % where each term is a kronecker product of $m$ matrices and $\mathbf{\Lambda}_{\mathcal{M}}$
        % occurs at the $j^{\mathrm{th}}$ location.
        %
        % We can show that the diagonal matrix $\mathbf{\Lambda}_{\mathcal{M},i}$ contains
        % $m^i$ non-zero elements  $(\mathbf{\Lambda}_{\mathcal{M},i})_{j(\mathbf{k})}
        % = \lambda_{k_1}+\dots+ \lambda_{k_i},\quad \mathbf{k}\in\mathbb{N}^i$ and ${j(\mathbf{k})}$
        % represents the lexicographical bijective indexing of $i$-tuples taking values
        % from $1,\dots,m$ and is given by |the combinator function.|
        disp(['Computing autonomous whisker at order ' num2str(i)])
        combinations = combinator(m,i,'p','r');
        Lambda_Mi = sum(Lambda_M(combinations),2);
        %%
        % Spectrum of $\mathcal{C}_i:=\left[\left(\mathbf{\Lambda}_{\mathcal{M},i}^{\top}\otimes\mathbf{B}\right)-\left(\mathbf{I}_{m^{i}}\otimes\mathbf{A}\right)\right]$
        %%
        % The matrix that needs to be inverted for solving the coefficients at the $i^{\mathrm{th}}$
        % order is given by
        %
        % $$\mathcal{C}_i:=\left[\left(\mathbf{\Lambda}_{\mathcal{M},i}^{\top}\otimes\mathbf{B}\right)-\left(\mathbf{I}_{m^{i}}\otimes\mathbf{A}\right)\right]$$
        %% Assemble RHS
        % $$\mathbf{L}_i =\sum_{j=2}^{l}\mathbf{F}_{j}\sum_{|\mathbf{p}|=i}\mathbf{S}_{p_{1}}\otimes\dots\otimes\mathbf{S}_{p_{j}}-\mathbf{B}\sum_{j=2}^{i-1}\mathbf{S}_{j}\sum_{|\mathbf{p}|=1}\mathbf{R}_{i+1-j}^{p_{1}}\otimes\dots\otimes\mathbf{R}_{i+1-j}^{p_{j}}$$
        %
        % where $l=\min(i,\Gamma_F)$ since we need to compute the summation at most
        % up to the order of the nonlinearity in $\mathbf{F}$.
        SIZE = [N, m*ones(1,i)];
        FS = sptensor(SIZE);
        %%
        % First term
        l = min(i,length(F));
        for j = 2:l           % Outer for loop can be parallelized - l cores
            % find values for j positive numbers summing up to i
            P = nsumk(j,i,'positive');
            FS = FS + tensor_composition(F{j},W_0,P,SIZE);
        end
        %%
        % Second term
        SR = sptensor(SIZE);
        % R_0i = R_0(i:-1:2); % used for parfor
        for j = 2:i-1         % Outer for loop can be parallelized - i-1 cores
            P = ones(j,j) + eye(j,j);
            R_j = {sptensor(speye(m,m)),R_0{i+1-j}};
            % R_j = {sptensor(speye(m,m)),R_0i{j}}; % used for parfor
            SR = SR + tensor_composition(W_0{j},R_j,P,SIZE);
        end
        L_i = FS - ttm(SR,B,1);
        %%
        % Convert $\mathbf{L}_i$ object to sparse vector
        L_i = sptenmat(permute(L_i,[1, ndims(L_i):-1:2]), 1:ndims(L_i));
        if isempty(L_i.vals)
            L_i = sparse(N_i,1);
        else
            L_i = sparse(L_i.subs(:,1),L_i.subs(:,2),L_i.vals,N_i,1);
        end
        %% *Solving for SSM coefficients and Reduced dynamics*
        % $${{\mathcal{C}}_{i}}~\mathfrak{vec}\left(\mathbf{S}_{i}\right)=\mathfrak{vec}\left(\mathbf{L}_{i}\right)-\underbrace{\left(\mathbf{I}_{m^{i}}\otimes\mathbf{B}\mathbf{S}_{1}\right)}_{{\mathcal{D}}_{i}}\mathfrak{vec}\left(\mathbf{R}_{i}\right)$$
        %
        W_0i = zeros(N,m^i); % W_0i is dense in nature
        R_0i = zeros(m,m^i); % This could be sparse
        
        L_i = reshape(L_i,N,[]);
        nRes = 0;
        paramStyle = obj.Options.paramStyle;
        parfor l = 1:m^i
            lambda_l = Lambda_Mi(l);
            C_l = lambda_l * B - A;
            L_il = L_i(:,l);
            %%
            % Checking for near-inner resonances
            J = find(abs(lambda_l - Lambda_M)<abstol);
            
            if ~isempty(J)                
                switch paramStyle
                    case 'normalform'
                        %%
                        % Choosing reduced dynamics using (near-)kernel of $\mathcal{C}_i$
                        R_0il = zeros(m,1); % for slicing use
                        for j = J
                            w_j = W_M(:,j);
                            % R_0i(j,l) = w_j'*L_il;
                            R_0il(j) = w_j'*L_il;
                        end
                        R_0i(:,l) = R_0il;
                    case 'graph'
                        R_0i(:,l) = W_M'*L_il;
                end                
                b_l = L_il - B * V_M * R_0i(:,l);
            else
                b_l = L_il;
            end
            nRes = nRes + numel(J);
            %%
            % Obtaining minimum-norm solution for $\mathbf{S}_i$ using $\texttt{lsqminnorm}$
            % which performs a complete orthogonal decomposition and is better suited for
            % sparse matrices as opposed to the Moore-Penrose pseudo-inverse ($\texttt{pinv}$).
            % We would like to use better iterative procedures moving forward, currently lsqlin
            % is not suited for complex data entries.
            W_0i(:,l) = lsqminnorm(C_l,b_l);
        end
        disp([num2str(nRes) ' (near) inner resonance(s) detected at order ' num2str(i)])
                
        W_0i = reshape(sptensor(W_0i(:)), [N, m*ones(1,i)]);
        R_0i = reshape(sptensor(R_0i(:)), [m, m*ones(1,i)]);
        W_0i = permute(W_0i,[1, ndims(W_0i):-1:2]);
        R_0i = permute(R_0i,[1, ndims(R_0i):-1:2]);
        
        multi_input =  [];
        
    case 'multiindex'
        
        k = i;           % convention here: call highest order k instead of i, call SSM dim l instead of m
        A   = obj.System.A; % A matrix
        B   = obj.System.B; % B matrix
        N   = obj.dimSystem; % Full system dimensionality in first-order form
        l   = numel(obj.E.spectrum); % dim(M): M is the master modal subspace
        F   = obj.System.F; % Full system Nonlinearity coefficients at different orders
        
        W_M = multi_input.W_M; % Right eigenvectors of the modal subspace
        H   = multi_input.H;    % composition coefficients
        %% Setup for case with symmetries
        % Conjugate center index at all orders
        
        z_k             = multi_input.Z_cci(k);     % Highest index that coefficients are computed for in conjugate ordering.
        % This is precisely the conjugate center inder (cci) at order k
        Lambda_M_vector = multi_input.Lambda_M_vector; % Vector with master evals sorted by size and imag/real
        K               = flip(sortrows(nsumk(l,k,'nonnegative')).',2);
        K               = K(:,multi_input.revlex2conj{k});
        K               = K(:,1:z_k);           % Set of order k multi-indices in conjugate ordering up to conjugate center index
        
        %% Make input into functions more clear
        % To make the code picture more clear we unify some input parameters into a
        % field.
        multi_input.N = N;
        multi_input.K = K;
        %% Assemble RHS
        % The right hand side terms can be split into three groups, which are be implemented
        % independently. Mixed Terms
        L_k = B * mixed_terms_in_SR(W_0,R_0,multi_input);
        
        % Force composition terms
        
        % The composition coefficients of power series
        H_k  = composition_coefficients_of_powerseries(W_0,H,multi_input);
        H{k} = H_k;
        % Now the force contribution for the equation at order k is computed
        L_k = L_k + autonomous_force_multi_index(F,H,multi_input);
        L_k = reshape(L_k,N*z_k,1);
        
        % Extract the near kernel of the coefficient matrix
        % coordinate directions do not change - we use the evals as in rev. lex
        % ordering - lambda_i has to be multiplied with i-th entry of a multi-index
        Lambda_Mk_vector         = sum(K.*Lambda_M_vector);
        [K_k,G_k,innerresonance] = explicit_kernel_multi_index(z_k,Lambda_M_vector,Lambda_Mk_vector,W_M,obj.Options.reltol);
        
        if innerresonance
            % here we use S_1 in rev_lex order since G_k is  constructed in rev_lex order for
            % correct reduced dynamics (coord directions)
            Skron = kron(speye(z_k),B*restore_full_coeff(W_0{1},[],multi_input.Z_cci(1),multi_input.conj2revlex{1}));
            R_0i   = G_k.' * K_k' * L_k;
            L_k   = L_k - Skron*R_0i;
        else
            R_0i   = sparse(l*z_k,1);
        end
        
        W_0i    = zeros(N,z_k);
        L_k    = reshape(L_k,N,z_k);
        
        % Solve the linear system for the SSM-coefficients
        parfor f = 1:z_k
            C_k        = B*Lambda_Mk_vector(f)-A;
            W_0i(:,f) = lsqminnorm(C_k,L_k(:,f));
        end
        R_0i        = reshape(R_0i,l,[]);
        H_k(:,:,1) = W_0i;
        
        %pass on composition coefficients
        H{k}       = H_k;
        multi_input.H = H;
end
% estime memory consumption from all variables in the current workspace
obj.solInfo.memoryEstimate(i) = monitor_memory('caller');
end