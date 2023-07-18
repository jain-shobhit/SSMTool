function [DX,DS] = compute_sensitivity_coefficients(obj,order,W,R,X)
% This function computes the derivatives of the reduced dynamics and SSM
% coefficients with respect to the forcing frequency.
% So far only implemented for periodic forcing.
%% System Properties

Omega  = obj.System.Omega;       % has to be column vector since kappas are stored in rows
A      = obj.System.A;           % A matrix
B      = obj.System.B;           % B matrix
N      = obj.dimSystem;          % full system size
W_M    = obj.E.adjointBasis ;    % Right eigenvectors of the modal subspace
l      = obj.dimManifold;        % dim(M): M is the master modal subspace
nKappa = obj.System.nKappa;

% Struct for passing variables to functions
field.N        = N;
field.l        = l;
field.ordering = 'revlex';
field.F_ord    = numel(obj.System.F);

% Struct for storing sensitivity coefficients
[DX,DS,kappas,field.Fext_ord] = struct_setup(obj,order);

% As these coefficients depend explicitly depend on omega
DX(1).Omega = Omega;
DS(1).Omega = Omega;

%% Zeroth order sensitivity equation

[W,R,field.H] = get_autonomous_coeffs(W,R); % reverse to lexicographical ordering
[X] = get_nonautonomous_coeffs(X); % reverse to lexicographical ordering


% coefficient matrix
for j = 1: nKappa
    
    % coefficient matrix
    C = A - 1i * kappas(j) * B;
    
    % right hand side
    RHS = (1i * kappas(j) * B) * X(j).W(1).coeffs;
    DX(j).W(1).coeffs = lsqminnorm(C,RHS);
    DX(j).W(1).ind = sparse(l,1);
end


%% Higher order sensitivity equation

parfor j = 1:nKappa
    
    kappa = kappas(j);
    for k = 1:order
        %Forcing and nonlinearity terms
        [G]  = Gnl(obj,field,k,DX(j));
        
        % Mixed Terms
        [XS] = RDX_plus_WDS(field,k,W,DX(j),R,DS(j));
        
        % Term from explicit omega derivative
        kappaX = 1i * kappa * X(j).W(k+1).coeffs;
        
        % Get resonant terms
        [E, I_k,K_lambda]   = resonant_terms(obj,k,kappa,'k');
        

        % Set reduced dynamics senitivity coefficients
        DS_jk               = sum( conj(W_M(:,E)).* ( G(:,I_k) - B*(XS(:,I_k))));
        DS(j).R(k+1).coeffs = sparse(E,I_k,DS_jk , l,nchoosek(k+l-1,l-1));
        DS(j).R(k+1).ind    = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
        
        %% Solve the sensitivity equation for the SSM sensitivity coefficients
    
        % Add DS order k contribution to the right hand side
        RHS                 = B* (XS +  coeffs_mixed_terms(k,1, W,DS(j).R,field,'R1') + kappaX) - G;
        
        for multi_idx = 1:nchoosek(k+l-1,l-1)
            C_ik = A - B * (K_lambda(multi_idx) + 1i * kappa*Omega); % Coefficient matrix
            DX(j).W(k+1).coeffs(:,multi_idx) = lsqminnorm(C_ik,RHS(:,multi_idx));
        end
        
        DX(j).W(k+1).ind = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
    end
    
    % Output coefficients in lexicographic ordering, with multi indices stored
    % in rows
    for k = 1:order+1 %index starts at 0
        DX(j).W(k) = coeffs_lex2revlex(DX(j).W(k),'TaylorCoeff');
        DS(j).R(k) = coeffs_lex2revlex(DS(j).R(k),'TaylorCoeff');
    end
end
end


function [W,R,H]               = get_autonomous_coeffs(W,R)
% Sets up the autonomous coefficients used in nonautonomous computation

%These quantities are all in lexicographic ordering, calculations are carried out in reverse
%lexicographic ordering. This is accounted for below.
W = coeffs_lex2revlex(W,'TaylorCoeff');
R = coeffs_lex2revlex(R,'TaylorCoeff');

%composition coefficients of power series
[H] = get_composition_coeffs(W);

end


function [X]               = get_nonautonomous_coeffs(X)
% Sets up the nonautonomous coefficients used in sensitivity computation

%These quantities are all in lexicographic ordering, calculations are carried out in reverse
%lexicographic ordering. This is accounted for below.
for i = 1:numel(X)
    X(i).W = coeffs_lex2revlex(X(i).W,'TaylorCoeff');
end
end

function [H]                     = get_composition_coeffs(W)
% This function reconstructs the composition coefficients for the computed
% SSM coefficients
%W_0 input in rev-lexicographic ordering, outputs H in rev-lexicographic ordering
field.ordering = 'revlex';

H = cell(1,numel(W));
H{1} = W(1).coeffs;

for k = 2:numel(W)
    field.k = k;
    H{k} = coeffs_composition(W,H,field);
end
end


function [G]                     = Gnl(obj,field,k,DX)
% Computes the forcing and nonlinearity contribution to the order k
% invariance equation for kappa_i

z_k = nchoosek(k+field.l-1,field.l-1);
G   = sparse(field.N,z_k);

for n = 2:k+1

    % NONLINEARITY
    % sum to k+1 since this term includes spatial derivatives
    if n <= field.F_ord && ~isempty(obj.System.F(n)) && ~isempty(obj.System.F(n).coeffs)
        G  = G + obj.System.F(n).coeffs* ...
            compute_sigma(obj.System.F(n).ind.',DX.W,k,field);
    end
    
end
end


function [XS]                    = RDX_plus_WDS(field,k,W,DX,R,DS)
% Computes the contributions of products of SSM and reduced dynamics
% coefficients to the order epsilon sensitivity equation

z_k = nchoosek(k+field.l-1,field.l-1);

RDX = sparse(field.N,z_k);
WDS = sparse(field.N,z_k);

% Terms with order 1 SSM coefficients (in epsilon)
field.mix = 'W1';
for m = 1:k %includes the zeroth order of W1
    if  ~isempty(DX.W(m).coeffs)
        RDX = RDX + coeffs_mixed_terms(k,m, DX.W, R,field,'W1');
    end
end

% Terms with order 1 reduced dynamics (in epsilon)
field.mix = 'R1';
for m = 2:k+1 % zeroth order in R1, no order k red. dyn.
    if ~isempty(DS.R(k-m+2).coeffs)
        WDS = WDS + coeffs_mixed_terms(k,m, W,DS.R,field,'R1');
    end
end

XS = WDS+RDX;
end


function [DX,DS,kappas,Fext_ord] = struct_setup(obj,order)
% Function that initialises the structs and some temporary arrays

l = obj.dimManifold;
N = obj.dimSystem;
nKappa = obj.System.nKappa;
k_kappa   = size(obj.System.Fext.data(1).kappa,1);
kappas    = zeros(k_kappa,nKappa);
% intitalise data structures to store coefficients
idle = repmat(struct('coeffs',[],'ind',[]),order+1  , 1);
DX  = repmat(struct('kappa' ,[],'W',idle,'Omega',[]),nKappa, 1);
DS  = repmat(struct('kappa' ,[],'R',idle,'Omega',[]),nKappa, 1);

Fext_ord = zeros(1,nKappa);

for i = 1:nKappa
    Fext_ord(i)  = numel(obj.System.Fext.data(i).F_n_k);
    kappa = obj.System.Fext.data(i).kappa;
    DX(i).kappa = kappa;
    DS(i).kappa = kappa;
    kappas(:,i)  = kappa;
    
    DX(i).W(1).coeffs = sparse(N,1);
    DX(i).W(1).ind    = sparse(l,1);
    DS(i).R(1).coeffs = sparse(l,1);
    DS(i).R(1).ind    = sparse(l,1);
end
end


function [E, I_k,K_lambda]       = resonant_terms(obj,k,kappa,order)
% This function finds the combinations of frequency multi-indices, master
% mode eigenvalues and the spatial multi-indices at zeroth and order k that
% lead to internal resonances.

Lambda = obj.E.spectrum;   % master modes eigenvalues
Omega  = obj.System.Omega; % forcing frequency
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
        
        K = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
        z_k = size(K,2);
        %vector with each element korresponding to summing multi_index k with all master lambdas
        K_lambda = sum(K .* Lambda);
        lambda_C_11 = repmat(Lambda,[1,z_k]) - repmat(K_lambda + 1i * (kappa*Omega),[l 1]);
        
        [E, I_k] = find(abs(lambda_C_11)<abstol); %I_k indicates the spatial multi-index the resonance corresponds to
end
end