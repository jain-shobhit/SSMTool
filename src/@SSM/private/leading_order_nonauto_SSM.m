function W1 = leading_order_nonauto_SSM(A,B,W_M,V_M,Lambda_M_vector,kappa_set,F_kappa,reltol,kNonauto,iNonauto,rNonauto,order,om)
% this function is adapted from compute_perturbed_wisker for the purpose of
% parallel computation. the reader may refer to the
% compute_perturbed_wisker for more details about this function.

% *Near external resonances and Leading order reduced dynamics*
F_kappacol = F_kappa(:);
[N,mm] = size(V_M);
numOm  = numel(om);
W1 = cell(numOm,1);
parfor j=1:numOm
    Omega = om(j);
    % *Near external resonances and Leading order reduced dynamics*
    [K, ~] = size(kappa_set);
    lambda_C_10 = 1i*repmat((kappa_set*Omega).',[mm 1]) - repmat(Lambda_M_vector,[1,K]);
    ref = min(abs(Lambda_M_vector));
    abstol = reltol * ref;
    % check external resonance
    [P, Q] = find(abs(lambda_C_10)<abstol);
    r_ext = length(Q);

    if r_ext    
        E_P = sparse(P, (1:r_ext).', true(r_ext,1), mm, r_ext );
        E_Q = sparse(Q, (1:r_ext).', true(r_ext,1), K, r_ext );
        W_P = W_M(:,P);
        G_10 = khatri_rao_product(E_Q,E_P)';
        K_10 = khatri_rao_product(E_Q,W_P);
        f_10 = G_10.' * (K_10' * F_kappacol);
        f_10 = reshape(f_10,mm,[]);
        l_10 = F_kappa - B * V_M * f_10;
    else
        l_10 = F_kappa;
        f_10 = [];
    end
    assert(~isempty(f_10), 'Near resonance does not occur, you may tune tol');
    f = f_10((kNonauto-1)*mm+2*iNonauto-1);
    assert(norm(f-rNonauto)<1e-3*norm(f), 'inner resonance assumption does not hold');

    % solving linear equations
    W1j = cell(1,order + 1);
    W10 = zeros(N,K);
    % Use conjugacy to reduce computations: fe^{i<kappa,omega>t} + \bar{f}e^{i<-kappa,omega>t}
    [redConj,mapConj] = conj_red(kappa_set, F_kappa);
    % Here redConj gives the reduced kappa sets by accounting for
    % complex conjugate, and mapConj is a cell array, where each entry maps to the
    % possible conjugacy pairs
    for jj = 1:numel(redConj)
        C_j  =  1i*dot(Omega,kappa_set(redConj(jj),:))*B - A;
        W10j = lsqminnorm(C_j,l_10(:,redConj(jj)));
        mapj = mapConj{jj};
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
    W1j{1}.coeffs = W10; W1j{1}.kappas = kappa_set;
    W1{j} = W1j;
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

