function [W1, R1] = compute_perturbed_whisker(obj, order,W0,R0,varargin)
% COMPUTE_PERTURBED_WHISKER This function computes the non-autonomous SSM   
%
% up to order 'order'.
%
% [W1, R1] = COMPUTE_PERTURBED_WHISKER(obj, order,W0,R0,varargin)
%
% obj:      SSM class object
% order:    approximation order up until which SSM is computed
% W0:       autonomous SSM coefficients
% R0:       autonomous RD coefficients
% varargin: if not stored in the obj class, this can be used to input Omega
%
% W1:       non-autonomous SSM coefficients
% R1:       non-autonomous RD coefficients
%
% See also: COMPUTE_WHISKER

%% Non-autonomous (quasi)periodic perturbation to whiskers of invariant manifolds
% We consider the mechanical system
%
% $$$\mathbf{B}\dot{\mathbf{x}} =\mathbf{Ax}+\mathbf{G}(\mathbf{x})+\epsilon\mathbf{F}(\mathbf{\phi},
% \mathbf{x}),$$$$\dot{\mathbf{\phi}}	=\mathbf{\Omega}$$
%
% with quasi-periodic forcing.
%
% In the non-autonomous setting, the SSM and the corresponding reduced dynamics
% would be parameterized by the angular variables $\mathbf{\phi}$, as well. In
% general, we may write
%
% $$    \mathbf{S(p,{\mathbf{\phi}})} = \mathbf{T}\mathbf{(p)} + \epsilon \mathbf{U}\mathbf{(p,\mathbf{\phi})}
% + O(\epsilon^2),$$
%
% $$    \mathbf{R(p,{\mathbf{\phi}})} = \mathbf{P}\mathbf{(p)} + \epsilon \mathbf{Q}\mathbf{(p,\mathbf{\phi})}
% + O(\epsilon^2),$$
%
% where $\mathbf{T}(\mathbf{p}),\mathbf{P}(\mathbf{p})$ recover the SSM and
% reduced dynamics coefficients in the unforced limit of $\epsilon=0.$
%
% These functions as well as the nonlinearity and the forcing are expanded in
% phase space coordinates. The time dependent coefficients of those expansions
% are furthermore expanded as a Fourier-series. As an example, the Force and the
% non-autonomous SSM-coefficients are given as
%
% $$    \mathbf {F}(\mathbf{x},\mathbf{\phi}) =     \left[  f^1(\mathbf{x},\mathbf{\phi}),
% \cdots ,   f^{2n}(\mathbf{x},\mathbf{\phi})     \right]^T,     \ f^i(\mathbf{x},\mathbf{\phi})
% = \sum_{\mathbf{n}\in \mathbb{N}^{2n}} F^i_{\mathbf{n}}(\mathbf{\phi}) \mathbf{x}^\mathbf{n}$$
%
% $$F^b_{\mathbf{k}}(\mathbf{\phi}) = \sum_{\mathbf{\eta} \in \mathbb{Z}^k}
% F^b_{\mathbf{k},\mathbf{\eta} } e^{i\langle \mathbf{\eta},\mathbf{\phi}\rangle}$$
%
% $$    \mathbf {U}(\mathbf{p},\mathbf{\phi}) =     \left[  u^1(\mathbf{p},\mathbf{\phi}),
% \cdots ,   u^{2n}(\mathbf{p},\mathbf{\phi})     \right]^T,     \ u^i(\mathbf{p},\mathbf{\phi})
% = \sum_{\mathbf{m}\in \mathbb{N}^{l}} U^i_{\mathbf{m}}(\mathbf{\phi}) \mathbf{p}^\mathbf{m}$$
%
% $$U^i_{\mathbf{k}}(\mathbf{\phi}) = \sum_{\mathbf{\eta} \in \mathbb{Z}^k}
% U^i_{\mathbf{k},\mathbf{\eta} } e^{i\langle \mathbf{\eta},\mathbf{\phi}\rangle}$$
%
% This leads to the invariane equation
%
% $$ \mathbf{B}  \bigg( \text{D}_\mathbf{p}( \mathbf{T}\mathbf{(p)})\mathbf{Q}\mathbf{(p,\mathbf{\phi})}
% + (\partial_\mathbf{p} \mathbf{U}\mathbf{(p,\mathbf{\phi)})  \mathbf{P}\mathbf{(p)}+
% (\partial_\mathbf{\phi} \mathbf{S(p,\mathbf{\phi})})   \mathbf{\Omega }  \bigg)=\\
% \ \ \ \ \mathbf{A}\mathbf{U}\mathbf{(p,\mathbf{\phi})} +     \big[\text{D}_\mathbf{x}\mathbf{G}
% \circ \mathbf{T}(\mathbf{p}) \big]\mathbf{U}(\mathbf{p},\mathbf{\phi})+
% \mathbf{F} (\mathbf{\phi},\mathbf{S(p,\mathbf{\phi})})$$
%
% The various expansions are plugged into this equation and then the equation
% is iteratively solved for the coefficients. The functions $\mathbf{T(p),P(p)}$
% are already known from the autonomous computation, their coefficients given
% by $\texttt{W0}$ (SSM) and $\texttt{R0}$ (reduced dynamics) .
%
% The external force $\mathbf{F}$ is input as a field of the property $\texttt{System}$
% of the SSM object. Since the equations for different frequency multi-indices
% decouple and the code is parallelised over these decoupled equations we want
% to make the read out hirarchy such that the first parameter corresponds to the
% frequency multi-indices. $\texttt{System.Fext.data(i)}$ indices into the $\texttt{i}$-th
% component of the field $\texttt{data}$ which is a struct array containing struct
% arrays. There is one such  contained struct array for each frequency multi-index.
%
% Every struct now contains two arrays with the coefficients and the spatial
% multi-indices respectively. The coefficients of the force for the $i$-th frequency-
% and ther order $k$ spatial multi-indices and their coefficients are stored in
% the rows of  $\texttt{data(i).F\_n\_k(k).ind} $ and the columns of $\texttt{data(i).F\_n\_k(k).coeffs}$
% . The $i$-th frequency multi-index is stored in $\texttt{data(i).kappa}$.
%
% The non-autonomous SSM and reduced dynamics coefficients are stored analogously.
% In $\texttt{W\_1(i).W(k).coeffs}$ the coefficients of the SSM expansion corresponding
% to $\mathbf{\eta}_i$ and order $k$ spatial multi-indices are stored. During
% the computation the multi-indices are stored in the columns of $\texttt{W\_1(i).W(k).ind}$
% in reverse lexicographic order, upon outputting the resulting coefficients however
% the storing scheme is reversed, in the output the multi-indices are stored in
% the rows in lexicographic ordering, the standard way of storing used throughout
% the software package.
%
% While in the documentation the frequency multi-indices are called $\mathbf{\eta}$
% for good distinguishability from spatial multi-indices in the code they are
% called $\texttt{kappa}$.
%% System Properties

if isempty(varargin) % FRC_po
    Omega = obj.System.Omega;       % has to be column vector since kappas are stored in rows
else % Parallelised FRC_lvlset
    Omega = varargin{1};
end

A      = obj.System.A;           % A matrix
B      = obj.System.B;           % B matrix
N      = obj.dimSystem;          % full system size
W_M    = obj.E.adjointBasis ;    % Right eigenvectors of the modal subspace
V_M    = obj.E.basis;            % Left eigenvectors of the modal subspace
l      = obj.dimManifold;        % dim(M): M is the master modal subspace
nKappa = obj.System.nKappa;
solver = obj.Options.solver;

% Struct for passing variables to functions
data.N        = N;
data.l        = l;
data.ordering = 'revlex';
data.F_ord    = numel(obj.System.F);


% Structs for storing coefficients
% each column in kappas corresponds to one kappa
if obj.Options.contribNonAuto % whether to ignore higher order
    [W1,R1,kappas,data.Fext_ord] = struct_setup(obj,order);
else %only zeroth order nonautonomous coefficients
    [W1,R1,kappas,data.Fext_ord] = struct_setup(obj,0);  
end 

% As these coefficients depend explicitly on omega at higher orders save it
% as property
W1(1).Omega = Omega;
R1(1).Omega = Omega;

%% Solving for coefficients with k=0
% The coefficient equation for this case reads
%
% $$\sum_{i=1}^{2n}     \underbrace{    \big(        \mathbf{(A)}_{bi}  - i
% \langle \mathbf{\eta},\mathbf{\Omega }\rangle        \mathbf{B}_{bi}    \big)
% }_{:= (\mathcal{L}_{\mathbf{0},\mathbf{\eta})_{bi}}}    U^i_{\mathbf{0},\mathbf{\eta}}=\sum_{j=1}^{l}
% (\mathbf{Bv}_j)_b      Q^j_{\mathbf{0},\mathbf{\eta}}-    F_{b,\mathbf{0},\mathbf{\eta}}$$

% Finding all force contributions at zeroth order
[F_0, idx_0] = zeroth_order_forcing(obj);
%% Find resonant terms

[ev_idx, harm_idx,~] = resonant_terms(obj,[],kappas(:,idx_0),Omega,'zero'); % contains harmonic idx.
r_ext    = length(ev_idx);
if obj.System.order == 1

%% Set reduced dynamics
% These sets now determine the bases of the near left kernel of the coefficient
% matrices $\mathbf{\mathcal{L}_{0,\eta}\cdot$ onto which we project the RHS of
% the coefficient equation and set it equal zero. This gives the explicit expression
% for the reduced dynamics coefficients.
%
% $$    Q^{e_i}_{\mathbf{0},\mathbf{\eta}_{f_i}}    =    \langle     \mathbf{w}_{e_i}
% ,     	\big[     		F^1_{\mathbf{0},\mathbf{\eta}_{f_i}}      		     		\cdots
% F^{2n}_{\mathbf{0},\mathbf{\eta}_{f_i}}     	\big]^T      \rangle$$

if r_ext
    Q_0 = sparse(ev_idx,harm_idx,sum(conj(W_M(:,ev_idx)).* F_0(:,idx_0(harm_idx))), l,numel(idx_0));
    RHS = B*V_M*Q_0 - F_0(:,idx_0);
else
    Q_0 = sparse(l,nKappa);
    RHS = - F_0(:,idx_0);
end
%% Solve for the order zeroth order SSM-coefficients W1, R1

run_idx = 1;
for j= idx_0    

    R1(j).R(1).coeffs = Q_0(:,run_idx);
    R1(j).R(1).ind    = sparse(l,1);    
    run_idx = run_idx + 1;
end

else % 2nd order computation
    M      = obj.System.M;
    C      = obj.System.C;
    THETA  = obj.E.adjointBasis(1:(N/2),:);
    PHI    = obj.E.basis(1:(N/2),:);
    Lambda_M_vector = obj.E.spectrum;
    V0     = sparse(N/2,1);
    
    if r_ext
        [harm_idx_un, ~,i_E_un] = unique(harm_idx.'); % unique harmonic positions
        
        Y0_all = sparse(N/2,numel(idx_0));
        ji = 1;
        for un_harm = harm_idx_un % Loop over all unique resonant harmonics
            un_ev = ev_idx((i_E_un == ji)); % All resonant eigenvalues for this harmonic
            
            Y0     = - F_0(1:(N/2),un_harm);
            Y0_all(:,un_harm) = Y0;
            Lambda_0_Om =  1i * (kappas(:,idx_0(un_harm))*Omega);
            
            [S0] = reduced_dynamics_second_order(un_ev,[],THETA,PHI,Lambda_M_vector,Lambda_0_Om, M,C,V0,Y0,l,1,0);
            
            R1(idx_0(un_harm)).R(1).coeffs = sparse(S0);
            R1(idx_0(un_harm)).R(1).ind    = sparse(l,1);
            ji = ji +1;
        end
    end
end


if obj.Options.contribNonAuto % whether to ignore higher order
    kappas_0 = kappas(:,idx_0);
    [redConj,mapConj] = conj_red(kappas_0, F_0(:,idx_0));
    
    % Make use of symmetry of leading order non-autonomous coefficients
    for j = 1:numel(redConj)
        
        if obj.System.order == 1
            C_j  =  A - 1i*dot(Omega,kappas_0(:,redConj(j)))*B ;
            W10j = solveinveq(C_j,RHS(:,redConj(j)),solver ); %lsqminnorm(C_j,RHS(:,redConj(j)));
            mapj = mapConj{j};
            
        else
            S0 = R1(j).R(1).coeffs;
            if isempty(S0)
                S0 = 0;
            end
            
            Lambda_0_Om =  1i * (kappas(:,redConj(j))*Omega);
            
            L_k = ( M * ((Lambda_0_Om + Lambda_M_vector.') .* PHI) + C*PHI ) * S0;
            L_k = L_k +  Y0_all(:,redConj(j));
            
            C_k = -(obj.System.K + Lambda_0_Om*C + Lambda_0_Om^2 *M );
            
            w_10    = solveinveq(C_k,L_k,solver);%lsqminnorm(C_k,L_k);

            
            w_10dot = Lambda_0_Om * w_10 + PHI * S0 + V0;
            
            W10j = [w_10;w_10dot];
            mapj = mapConj{j};
        end
        
        switch numel(mapj)
            case 1
                W1(idx_0(mapj(1))).W(1).coeffs = W10j;
                W1(idx_0(mapj(1))).W(1).ind    = sparse(l,1);
            case 2

                W1(idx_0(mapj(1))).W(1).coeffs = W10j;
                W1(idx_0(mapj(2))).W(1).coeffs = conj(W10j);
                
                W1(idx_0(mapj(1))).W(1).ind    = sparse(l,1);
                W1(idx_0(mapj(2))).W(1).ind    = sparse(l,1);
            otherwise
                error('there exist redundancy in kappa of external forcing');
        end
    end
    %% Solving for coefficients with k>0
    % The coefficient equation for this case reads
    %
    % $$\sum_{i=1}^{2n}     \underbrace{     \bigg(       \mathbf{(A)}_{bi}
    % -        \mathbf{B}_{bi}         \big[            \sum_{j=1}^l k_j  \lambda_{j}
    % +            i\langle \mathbf{\Omega}, \mathbf{\eta} \rangle         \big]
    % \bigg)     }_{:= (\mathcal{L}_{\mathbf{k},\mathbf{\eta}})_{bi}}    U^i_{\mathbf{k},\mathbf{\eta}}\\=
    % \sum_{i=1}^{2n} \mathbf{B}_{bi}\sum_{j=1}^l    \bigg[          \sum_{\mathbf{m},
    % \mathbf{u}\in \mathbb{N}^l , \ \mathbf{m+u} - \mathbf{\hat{e}}_j = \mathbf{k}}
    % m_j        T^i_{\mathbf{m}}        Q^j_{\mathbf{u},\mathbf{\eta}}        +
    % \sum_{\mathbf{m,u} \in \mathbb{N}^l, \ |\mathbf{m}|<k \ \ \mathbf{m+u} - \hat{\mathbf{e}}_j
    % = \mathbf{k}} m_j U^i_{\mathbf{m},\mathbf{\eta}} P^j_{\mathbf{u}}    \bigg]
    % \\  \ \ \ -    \sum_{\mathbf{n}\in \mathbb{N}^{2n}, |\mathbf{n}|<k}
    % F^b_{\mathbf{n},\mathbf{\eta}} \pi_{\mathbf{n,k}}-    \sum_{\mathbf{n}\in \mathbb{N}
    % ^{2n}, \ \            |\mathbf{n}| \geq 2}       G^b_{\mathbf{n}}\sigma_{\mathbf{k},
    % \mathbf{n}, \mathbf{\eta}}$$
    
    % Get autonomous coefficients and composition coefficients in rev. lex.
    % ordering
    if order>0
        [W0,R0,data.H] = get_autonomous_coeffs(W0,R0);
    end
    %% Perform Nonautonomous Calculation
    % We loop over all orders of spatial multi-indices. Within that there is a loop
    % over all the frequency multi-indices.
    
    for i = 1:nKappa
        
        for k = 1:order
            %% Calculating the RHS
            %
            
            %Forcing and nonlinearity terms
            [FG]  = Fext_plus_Gnl(obj,data,i,k,W1(i));
            
            % Mixed Terms
            [WR] = W1R0_plus_W0R1(data,k,W0,W1(i),R0,R1(i));
            %% Find resonant terms
            
            [ev_idx, multi_idx,Lambda_K]   = resonant_terms(obj,k,kappas(:,i),Omega,'k'); % F contains multi-index pos.
            %% Set reduced dynamics
            
            if obj.System.order == 1
                [R1_ik,W1_ik] = first_order_computation(kappas,Omega,W_M,FG,WR,A,B,i,k,ev_idx,multi_idx,l,data,Lambda_K,W0,solver);  
            else
                [R1_ik,W1_ik] = second_order_computation(kappas,Omega,FG,WR,i,k,ev_idx,multi_idx,data,M,C,THETA,PHI,Lambda_K,Lambda_M_vector,obj.System.K,solver);
            end
            
            R1(i).R(k+1) = R1_ik(k+1);
            W1(i).W(k+1).coeffs = W1_ik.coeffs;
            
            if l >1
                W1(i).W(k+1).ind = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
            else
                W1(i).W(k+1).ind = k;
            end
        end
        
        % Output coefficients in lexicographic ordering, with multi indices stored
        % in rows
        for k = 1:order+1 %index starts at 0
            W1(i).W(k) = coeffs_lex2revlex(W1(i).W(k),'TaylorCoeff');
            R1(i).R(k) = coeffs_lex2revlex(R1(i).R(k),'TaylorCoeff');
        end
    end
    

end



end

function [W] = solveinveq(C,R,solver)
% solves C*W = R for W with given solver
switch solver
    case 'lsqminnorm'
        W = lsqminnorm(C,R);
    case 'linsolve'
        W = linsolve(C,R);
    case 'backslash'
        W = C \R;
    case 'pinv'
        W = pinv(C)*R;
    case 'inv'
        W = inv(C)*R;
end
end

function [R1_ik,W1_ik] = second_order_computation(kappas,Omega,FG,WR,i,k,I,F,data,M,C,THETA,PHI,Lambda_K,Lambda_M_vector,Ksys,solver)
l = data.l;
N = data.N;

Ym     = ( C * WR(1:(N/2),:) + M * WR((N/2+1):end,:)) - FG(1:(N/2),:);
Vm     = WR(1:(N/2),:);

Lambda_K_Om = Lambda_K + 1i * (kappas(:,i)*Omega);

z_k = nchoosek(k+l-1,l-1);

[Rk] = reduced_dynamics_second_order(I,F,THETA,PHI,Lambda_M_vector,Lambda_K_Om, M,C,Vm,Ym,l,z_k,k);


w_1i    = zeros(N/2,z_k);
w_1idot = zeros(N/2,z_k);

for f = 1:z_k
    
    L_k = ( M * ((Lambda_K_Om(f) + Lambda_M_vector.') .* PHI) + C*PHI ) * Rk(:,f);
    L_k = L_k + Lambda_K_Om(f)*M*Vm + Ym;
    
    C_k = -(Ksys + Lambda_K_Om(f)*C + Lambda_K_Om(f)^2 *M );
    

    w_1i(:,f)    = solveinveq(C_k,L_k(:,f),solver);

    w_1idot(:,f) = Lambda_K_Om(f) * w_1i(:,f) + PHI * Rk(:,f) + Vm(:,f);
end

W1_ik.coeffs       = [w_1i;w_1idot];
R1_ik(k+1).coeffs       = sparse(Rk);
if l>1
    R1_ik(k+1).ind    = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
else
    R1_ik(k+1).ind  = k;
end
end

function [R1_ik,W1_ik] = first_order_computation(kappas,Omega,W_M,FG,WR,A,B,i,k,I,F,l,data,K_lambda,W0,solver)
% Computes the SSM coeffs and reduced dynamics coefficients using first
% order implementation

R1_ik_coeff     = sum( conj(W_M(:,I)).* ( FG(:,F) - B*(WR(:,F))));
R1_ik(k+1).coeffs    = sparse(I,F,R1_ik_coeff , l,nchoosek(k+l-1,l-1));
if l > 1
    R1_ik(k+1).ind    = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
else
    R1_ik(k+1).ind = k;
end
%% Solve the coefficient equation for the SSM coefficients
% Add R1 order k contribution to the right hand side

RHS                 = B* (WR +  coeffs_mixed_terms(k,1, W0,R1_ik,data,'R1')) - FG;

for j = 1:nchoosek(k+l-1,l-1)
    C_i = A - B * (K_lambda(j) + 1i * kappas(:,i)*Omega); % Coefficient matrix
    W1_ik.coeffs(:,j) = solveinveq(C_i,RHS(:,j),solver);
end
end

function [redConj,mapConj] = conj_red(kappa_set,F_kappa)
% This function detects complex conjugate relations between forcing. For instance,
% when kappa_set = [1,-1,2,3,-3] and F_kappa = [1;1;2;3;4], it will return
% redConj = [1,3,4,5] with mapConj = {[1 2],3,4,5}
redConj = [];
mapConj = [];
assert(numel(kappa_set)==numel(unique(kappa_set)),'there exist redundancy in kappa of external forcing');
kappa = kappa_set;
while ~isempty(kappa)
    ka = kappa(1);
    ka_redConj = find(kappa_set==ka);
    redConj = [redConj;ka_redConj];
    % find the conjugate one if it exists
    ka_conj = find(kappa_set==-ka);
    if ~isempty(ka_conj) && norm(conj(F_kappa(:,ka_redConj))-F_kappa(:,ka_conj))<1e-6*norm(F_kappa(:,ka_conj))
        mapConj = [mapConj, {[ka_redConj,ka_conj]}];
        kappa = setdiff(kappa,[ka,-ka],'stable');
    else
        mapConj = [mapConj, {ka_redConj}];
        kappa = setdiff(kappa,ka,'stable');
    end
end
end
%%
%
%
%

function [W1,R1,kappas,Fext_ord] = struct_setup(obj,order)
% Function that initialises the structs and some temporary arrays

l = obj.dimManifold;
N = obj.dimSystem;
nKappa = obj.System.nKappa;
k_kappa   = size(obj.System.Fext.data(1).kappa,1);
kappas    = zeros(k_kappa,nKappa);
% intitalise data structures to store coefficients
idle = repmat(struct('coeffs',[],'ind',[]),order+1  , 1);
W1  = repmat(struct('kappa' ,[],'W',idle,'Omega',[]),nKappa, 1);
R1  = repmat(struct('kappa' ,[],'R',idle,'Omega',[]),nKappa, 1);

Fext_ord = zeros(1,nKappa);

for i = 1:nKappa
    Fext_ord(i)  = numel(obj.System.Fext.data(i).F_n_k);
    kappa = obj.System.Fext.data(i).kappa;
    W1(i).kappa = kappa;
    R1(i).kappa = kappa;
    kappas(:,i)  = kappa;
    
    W1(i).W(1).coeffs = sparse(N,1);
    W1(i).W(1).ind    = sparse(l,1);
    R1(i).R(1).coeffs = sparse(l,1);
    R1(i).R(1).ind    = sparse(l,1);
end
end
%%
%
%
%

function [F_0, idx_0]            = zeroth_order_forcing(obj)
% Finding all force contributions at zeroth order

nKappa = obj.System.nKappa;
N         = obj.dimSystem;
F_0 = zeros(N,nKappa);
for i = 1:nKappa
    if ~isempty(obj.System.Fext.data(i).F_n_k(1).coeffs)
        F_0(:,i)= obj.System.Fext.data(i).F_n_k(1).coeffs;  % each column corresponds to one kappa
    end
end
idx_0 = find(any(F_0~=0)); % index for all kappas that contribute at zeroth order
end

%%
%
%
%

function [W0,R0,H]               = get_autonomous_coeffs(W0,R0)
% Sets up the autonomous coefficients used in nonautonomous computation

%These quantities are all in lexicographic ordering, calculations are carried out in reverse
%lexicographic ordering. This is accounted for below.

W0 = coeffs_lex2revlex(W0,'TaylorCoeff');
R0 = coeffs_lex2revlex(R0,'TaylorCoeff');

%composition coefficients of power series
[H] = get_composition_coeffs(W0);

end
%%
%
%
%

function [H]                     = get_composition_coeffs(W0)
% This function reconstructs the composition coefficients for the computed
% SSM coefficients
%W_0 input in rev-lexicographic ordering, outputs H in rev-lexicographic ordering
field.ordering = 'revlex';

H = cell(1,numel(W0));
H{1} = W0(1).coeffs;

for k = 2:numel(W0)
    field.k = k;
    H{k} = coeffs_composition(W0,H,field);
end
end
%%
%
%
%

function [F]                     = Fext_plus_Gnl(obj,data,i,k,W1)
% Computes the forcing and nonlinearity contribution to the order k
% invariance equation for kappa_i

z_k = nchoosek(k+data.l-1,data.l-1);
Force  = sparse(data.N,z_k);
G_nl   = sparse(data.N,z_k);
if data.l > 1
    K   = flip(sortrows(nsumk(data.l,k,'nonnegative')).',2); %order k multi-indices
else
    K = k;
end

for n = 2:k+1
    % FORCING
    %sum to k+1 since index starts at 0 for k=0
    if  n <= data.Fext_ord(i) && ~isempty(obj.System.Fext.data(i).F_n_k(n).coeffs)
        
        F_coeff = obj.System.Fext.data(i).F_n_k(n).coeffs;
        F_ind   = obj.System.Fext.data(i).F_n_k(n).ind.';
        Force   = Force + F_coeff * compute_pi(F_ind,K, data);
        
    end
    % NONLINEARITY
    % sum to k+1 since this term includes spatial derivatives
    if n <= data.F_ord && ~isempty(obj.System.F(n)) && ~isempty(obj.System.F(n).coeffs)
        G_nl  = G_nl + obj.System.F(n).coeffs* ...
            compute_sigma(obj.System.F(n).ind.',W1.W,k,data);
    end
    
end
F = Force + G_nl;
end
%%
%
%
%

function [WR]                    = W1R0_plus_W0R1(data,k,W0,W1,R0,R1)
% Computes the contributions of products of SSM and reduced dynamics
% coefficients to the order epsilon invariance equation

z_k = nchoosek(k+data.l-1,data.l-1);

W1R0 = sparse(data.N,z_k);
W0R1 = sparse(data.N,z_k);

% Terms with order 1 SSM coefficients (in epsilon)
data.mix = 'W1';
for m = 1:k %includes the zeroth order of W1
    if  ~isempty(W1.W(m).coeffs)
        W1R0 = W1R0 + coeffs_mixed_terms(k,m, W1.W, R0,data,'W1');
    end
end

% Terms with order 1 reduced dynamics (in epsilon)
data.mix = 'R1';
for m = 2:k+1 % zeroth order in R1, no order k red. dyn.
    if ~isempty(R1.R(k-m+2).coeffs)
        W0R1 = W0R1 + coeffs_mixed_terms(k,m, W0,R1.R,data,'R1');
    end
end

WR = W0R1+W1R0;
end
%%
%
%
%

function [E, I_k,K_lambda]       = resonant_terms(obj,k,kappa,Omega,order)
% This function finds the combinations of frequency multi-indices, master
% mode eigenvalues and the spatial multi-indices at zeroth and order k that
% lead to internal resonances.

Lambda = obj.E.spectrum;   % master modes eigenvalues
l      = obj.dimManifold;
% Tolerance for resonances
ref = min(abs(Lambda));
abstol = obj.Options.reltol * ref;

switch order
    case 'zero'
        %% Find zeroth order resonant terms
        % We determine the near inner resonances of the coefficient matrix where
        %
        % $$     \lambda_{j} - i\langle \mathbf{\eta}_{f}, \mathbf{\Omega } \rangle\approx
        % 0, \ e\in \{1,...,l\}, \ f \in\{1,...,K\}$$
        %
        % holds. The index pairs that fulfill this condition are stored.
        %
        % $$E := \{ e_1,  ... ,e_{r_{ext}} \in \{1,...,l\}\} \\ F := \{ f_1,  ... ,f_{r_{ext}}
        % \in \{1,...,K\}\}$$
        
        % kappa in this case contains all kappas
        lambda_C_10 =  repmat(Lambda,[1,size(kappa,2)]) - 1i*repmat(kappa*Omega,[l 1]);
        
        [E, I_k] = find(abs(lambda_C_10)<abstol);
        K_lambda = [];
        
        % I_k contains the frequency index
    case 'k'
        %% Find higher order resonant terms
        % The coefficient matrix for frequency multi-index $\mathbf{\eta}$ shows singularities
        % if the resonance condition
        %
        % $$    \lambda_e - \bigg( \sum_{j=1}^l k_j\lambda_j     + i \langle \mathbf{\Omega},
        % \mathbf{\eta} \rangle \bigg) \approx 0$$
        %
        % is fulfilled for some $\lambda_e$ in the master subspace. We therefore have
        % to find all such resonant combinations.
        
        %Find the resonances
        
        if l > 1
            K = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
        else
            K = k;
        end
        z_k = size(K,2);
        %vector with each element korresponding to summing multi_index k with all master lambdas
        K_lambda = sum(K .* Lambda);
        lambda_C_11 = repmat(Lambda,[1,z_k]) - repmat(K_lambda + 1i * (kappa*Omega),[l 1]);
        
        [E, I_k] = find(abs(lambda_C_11)<abstol); %I_k indicates the spatial multi-index the resonance corresponds to
end
end


function [Rk] = reduced_dynamics_second_order(I,F,THETA,PHI,Lambda,Lambda_K, M,C,Vm,Ym,l,z_k,order)
% I         - Eigenvalue positions that are resonant
% F         - Multi-index positions that are resonant for order > 0,
%               otherwise the frequency harmonic
% Lambda    - vector containing master spectrum
% Lambda_K  - sum(Lambda * K + kappa * Omega)
% M         - Mass matrix
% C         - Damping Matrix
% Vm        - Velocity part of the lower order and forcing RHS
% Ym        - Displacement part of the lower order and forcing RHS
% l         - SSM dimension
% z_k       - number of multi-indices at order k
% order     - current order of computation


Rk = zeros(l,z_k);

if order == 0
    
    THETA_I = THETA(:,I);
    
    RHS   = -THETA_I' * (Lambda_K*M*Vm +Ym);
    C_0_r = eye(size(RHS,1)); % Fix leading order coefficient to exact resonance    
    Rk(I,:) = lsqminnorm(C_0_r,RHS);
    
else
    % unique multi indices
    if any(F)
        [F_un, ~,i_F_un] = unique(F.');
        
        % Loop over multi-indices that lead to resonance
        fi = 1;
        
        for f = F_un
            I_f = I((i_F_un == fi)); % All resonant eigenvalues for this multi - index
            THETA_f = THETA(:,I_f);
            PHI_f   = PHI(:,I_f);

            RHS   = - Lambda_K(f).*THETA_f' *M*Vm(:,f) -THETA_f'*Ym(:,f);
            %C_k_r =  THETA_f' * ( C* PHI_f + M * ((Lambda_K(f) + Lambda(I_f).') .* PHI_f)); %Coefficient Matrix
            C_k_r = eye(size(RHS,1));
            
            Rk(I_f,f) = lsqminnorm(C_k_r,RHS);
         
            fi = fi +1;
        end
    end
    
end
end

