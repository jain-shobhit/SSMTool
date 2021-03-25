function [W1, R1] = compute_perturbed_whisker(obj, order)
%% COMPUTE_PERTURBED_WHISKER Non-autonomous quasiperiodic perturbation to whiskers of invariant manifolds
% We currently perform only leading order calculations $$\mathbf{B}\dot{\mathbf{z}} 
% =\mathbf{Az}+\mathbf{F}(\mathbf{z})+\epsilon\mathbf{F}^{ext}(\mathbf{\phi}),$$
% 
% $$\dot{\mathbf{\phi}}	=\mathbf{\Omega}$$
% 
% In the non-autonomous setting, the SSM and the corresponding reduced dynamics 
% would be parameterized by the angular variables $\mathbf{\phi}$, as well. In 
% general, we may write
% 
% $$\mathbf{S}(\mathbf{p},\mathbf{\phi})=\mathbf{S}_{0}(\mathbf{p})+\epsilon\mathbf{S}_{1}(\mathbf{p},\mathbf{\phi})+\mathcal{O}\left(\epsilon^{2}\right),$$
% 
% $$\mathbf{R}(\mathbf{p},\mathbf{\phi})=\mathbf{R}_{0}(\mathbf{p})+\epsilon\mathbf{R}_{1}(\mathbf{p},\mathbf{\phi})+\mathcal{O}\left(\epsilon^{2}\right),$$
% 
% where $\mathbf{S}_{0}(\mathbf{p}),\mathbf{R}_{0}(\mathbf{p})$ recover the 
% SSM and reduced dynamics coefficients in the unforced limit of $\epsilon=0.$ 
% The non-autonomous invariance equation is given as
% 
% $$\mathbf{B}\left[D\mathbf{S}_{0}(\mathbf{p})\mathbf{R}_{1}(\mathbf{p},\mathbf{\phi})+\partial_{\mathbf{p}}\mathbf{S}_{1}(\mathbf{p},\mathbf{\phi})\mathbf{R}_{0}(\mathbf{p})+\partial_{\mathbf{\mathbf{\phi}}}\mathbf{S}_{1}(\mathbf{p},\mathbf{\phi})\cdot\mathbf{\Omega}\right]=\left[\mathbf{A}+D\mathbf{F}(\mathbf{S}_{0}(\mathbf{p}))\right]\mathbf{S}_{1}(\mathbf{p},\mathbf{\phi})+\mathbf{F}^{ext}(\mathbf{\phi}).$$
% 
% At the leading order in $$\mathbf{p}$$(namely $$\mathbf{p}^0$$), this invariance 
% equation is given by
% 
% $$\mathbf{B}\left[\mathbf{S}_{0,1}\mathbf{R}_{1,0}(\mathbf{\phi})+\partial_{\mathbf{\mathbf{\phi}}}\mathbf{S}_{1,0}(\mathbf{\phi})\cdot\mathbf{\Omega}\right]=\mathbf{A}\mathbf{S}_{1,0}(\mathbf{\phi})+\mathbf{F}^{ext}(\mathbf{\phi}),$$
% 
% which is system of linear differential equations for the unknown, time-dependent 
% coefficients $$\mathbf{S}_{1,0}(\mathbf{\phi})$$, where we may choose the reduced 
% dynamics $$\mathbf{R}_{1,0}(\mathbf{\phi})$$ appropriately. Consider Fourier 
% expansion as follows
% 
% $$\mathbf{F}^{ext}(\mathbf{\Omega}t)=\sum_{\mathbf{\kappa}\in\mathbb{Z}^{k}}\mathbf{F}_{\mathbf{\kappa}}e^{\iota\left\langle\mathbf{\kappa},\mathbf{\Omega}\right\ranglet}$$
% 
% $$\mathbf{S}_{1,0}(\mathbf{\Omega}t)=\sum_{\mathbf{\kappa}\in\mathbb{Z}^{k}}\mathbf{s}_{1,0,\mathbf{\kappa}}e^{\iota\left\langle\mathbf{\kappa},\mathbf{\Omega}\right\ranglet}$$
% 
% $$\mathbf{R}_{1,0}(\mathbf{\Omega}t)=\sum_{\mathbf{\kappa}\in\mathbb{Z}^{k}}\mathbf{r}_{1,0,\mathbf{\kappa}}e^{\iota\left\langle\mathbf{\kappa},\mathbf{\Omega}\right\ranglet}$$
% 
% Upon comparing Fourier coefficients at order $\mathbf{\kappa}$, we get
% 
% $$\underbrace{[\iota\left\langle \mathbf{\kappa},\mathbf{\Omega}\right\rangle\mathbf{B}-\mathbf{A}]}_{\mathbf{\mathcal{A}}_{0,\mathbf{\kappa}}}\mathbf{s}_{1,0,\mathbf{\kappa}}=\mathbf{F}_{\mathbf{\kappa}}-\mathbf{B}\mathbf{S}_{0,1}\mathbf{r}_{1,0,\mathbf{\kappa}}.$$
% 
% Vectorized version is written as
% 
% $\mathbf{\mathfrak{C}}_{1,0}\mathbf{\mathfrak{s}}_{1,0}=\mathbf{\mathfrak{f}}_{1,0}-\mathbf{\mathfrak{D}}_{1,0}\mathbf{\mathfrak{r}}_{1,0}=:\mathbf{\mathfrak{l}}_{10}$,
% 
% where
% 
% $\mathbf{\mathfrak{C}}_{1,0}=\mathrm{diag}\left(\mathbf{\mathcal{A}}_{0,\mathbf{\kappa}_{1}},\dots,\mathbf{\mathcal{A}}_{0,\mathbf{\kappa}_{K}}\right)=\mathfrak{\mathbf{R}}_{0}\otimes\mathbf{B}-\mathbf{I}_{K}\otimes\mathbf{A}\in\mathbb{C}^{KN\times 
% KN}$,
% 
% $\mathfrak{R}_{1,0}:=\mathrm{diag}\left(\iota\left\langle \mathbf{\kappa}_{1},\mathbf{\Omega}\right\rangle 
% ,\dots,\iota\left\langle \mathbf{\kappa}_{K},\mathbf{\Omega}\right\rangle \right)\in\mathbb{C}^{KN\times 
% KN}$,
% 
% $\mathbf{\mathfrak{s}}_{1,0}:=\mathfrak{vec}([\mathbf{s}_{1,0,\mathbf{\kappa}_{1}},\dots,\mathbf{s}_{1,0,\mathbf{\kappa}_{K}}])\in\mathbb{C}^{KN}$,
% 
% $\mathbf{\mathfrak{D}}_{1,0}:=\mathbf{I}_{K}\otimes\mathbf{B}\mathbf{S}_{0,1}\in\mathbb{C}^{KN\times 
% KM}$,
% 
% $\mathbf{\mathfrak{r}}_{1,0}:=\mathfrak{vec}([\mathbf{r}_{1,0,\mathbf{\kappa}_{1}},\dots,\mathbf{r}_{1,0,\mathbf{\kappa}_{K}}])\in\mathbb{C}^{KM}$,
% 
% $$\mathbf{\mathfrak{f}}_{1,0}:=\mathfrak{vec}\left(\left[\mathbf{F}_{\mathbf{\kappa}_{1}},\dots,\mathbf{F}_{\mathbf{\kappa}_{K}}\right]\right).$$
% 
% This function returns the leading order SSM coefficients $$\mathbf{\mathfrak{s}}_{1,0}$$ 
% and reduced dynamics coefficients $\mathbf{\mathfrak{r}}_{1,0}$.
R1 = cell(1,order + 1);
kappa_set= obj.System.Fext.kappas; % each row corresponds to one kappa
F_kappa = obj.System.Fext.coeffs;  % each column corresponds to one kappa
Omega = obj.System.Omega;
A = obj.System.A;           % A matrix
B = obj.System.B;           % B matrix
N = obj.dimSystem;          % full system size
W_M = obj.E.adjointBasis;   % Right eigenvectors of the modal subspace
V_M = obj.E.basis;          % Left eigenvectors of the modal subspace
m = obj.dimManifold; % dim(M): M is the master modal subspace
% *Near external resonances and Leading order reduced dynamics*
Lambda_M_vector = obj.E.spectrum;
[K, k] = size(kappa_set);
lambda_C_10 = 1i*repmat((kappa_set*Omega).',[m 1]) - repmat(Lambda_M_vector,[1,K]);
ref = min(abs(Lambda_M_vector));
abstol = obj.Options.reltol * ref;
%% 
% For a system with $r_{ext}$ near-resonances with the external forcing, let 
% the index sets
% 
% $$\mathcal{P}	:=\{p_{1},\dots,p_{r_{ext}}\in\{1,\dots,M\}\},\mathcal{Q}	:=\{q_{1},\dots 
% q_{r_{ext}}\in\{1,\dots,K\}\}$$
% 
% be defined such that 
% 
% $$\left(\lambda_{p_{j}}-\iota\left\langle \mathbf{\kappa}_{q_{j}},\mathbf{\Omega}\right\rangle 
% \right)\approx 0\,\quad\forall\,j\in\{1,\dots,r_{ext}\}.$$
[P, Q] = find(abs(lambda_C_10)<abstol);
r_ext = length(Q);
%% 
% Getting reduced dynamics coefficient
% 
% The left near-kernel of $$\mathbf{\mathfrak{C}}_{0}$ is useful for normal 
% form parametrization of reduced dynamics
% 
% $$\mathbf{\mathfrak{K}}_{1,0}=\mathbf{E}_{\mathcal{Q}}\odot\mathbf{W}_{\mathcal{P}}\in\mathbb{C}^{NK\times 
% r_{ext}},$$
% 
% where
% 
% $$\mathbf{E}_{\mathcal{Q}}=\left[\mathbf{e}_{q_{1}},\dots,\mathbf{e}_{q_{r_{ext}}}\right]\in\mathbb{R}^{K\times   
% r_{ext}},\mathbf{W}_{\mathcal{P}}=\left[\mathbf{w}_{p_{1}},\dots,\mathbf{w}_{p_{r_{ext}}}\right]\in\mathbb{C}^{N\times   
% r_{ext}}$$.
% 
% The reduced dynamics is then given as 
% 
% $$\mathbf{\mathfrak{r}}_{1,0}=\mathbf{\mathfrak{G}}_{1,0}^{\top}\mathbf{\mathfrak{K}}_{1,0}^{\star}\mathbf{\mathfrak{f}}_{1,0}$,
% 
% where
% 
% $$\mathbf{\mathfrak{G}}_{1,0} := \left(\mathbf{E}_{\mathcal{Q}}\odot\mathbf{E}_{\mathcal{P}}\right)^{\top}\in\mathbb{C}^{r_{ext}\times   
% KM}$.
if r_ext    
    E_P = sparse(P, (1:r_ext).', true(r_ext,1), m, r_ext );
    E_Q = sparse(Q, (1:r_ext).', true(r_ext,1), K, r_ext );
    W_P = W_M(:,P);
    G_10 = khatri_rao_product(E_Q,E_P)';
    K_10 = khatri_rao_product(E_Q,W_P);
    f_10 = G_10.' * (K_10' * F_kappa(:));
    f_10 = reshape(f_10,m,[]);
    l_10 = F_kappa - B * V_M * f_10;
else
    l_10 = F_kappa;
    f_10 = [];
end
R1{1}.coeffs = f_10; R1{1}.kappas = kappa_set;

%% 
% Assembling Coefficient matrix and solving system.
if obj.Options.contribNonAuto % whether to ignore higher order
    W1 = cell(1,order + 1);
    W10 = zeros(N,K);
    % Use conjugacy to reduce computations: fe^{i<kappa,omega>t} + \bar{f}e^{i<-kappa,omega>t}
    [redConj,mapConj] = conj_red(kappa_set, F_kappa);
    % Here redConj gives the reduced kappa sets by accounting for
    % complex conjugate, and mapConj is a cell array, where each entry maps to the
    % possible conjugacy pairs
    for j = 1:numel(redConj)
        C_j  =  1i*dot(Omega,kappa_set(redConj(j),:))*B - A;
        W10j = lsqminnorm(C_j,l_10(:,redConj(j)));
        mapj = mapConj{j};
        switch numel(mapj)
            case 1
                W10(:,mapj) = W10j;
            case 2
                W10(:,mapj(1)) = W10j;
                W10(:,mapj(2)) = conj(W10j);
            otherwise
                error('there exist redundancy in kappa of external forcing');
        end
    end
    W1{1}.coeffs = W10; W1{1}.kappas = kappa_set;
else
    W1 = [];
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
