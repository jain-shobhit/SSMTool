function extract_FRS(obj,oid,modes,order,mFreqs,parRange,outdof,optdof,scale_state,scale_obs,varargin)
% EXTRACT_FRS This function extracts the frequency response surface of
% mechanical systems subject to periodic forcing via SSM-based ROMs. We
% perform two-dimensional continuation of fixed points of the ROMs to cover
% the frequency response surface.
%
% FRC = EXTRACT_FRS(OBJ,OID,RESONANT_MODES,ORDER,MFREQS,PARRANGE,OUTDOF,OPTDOF,VARARGIN)
%
% oid:      runid of continuation
% modes:    master subspace
% order:    expansion order of SSM
% mFreqs:   internal resonance relation vector
% parRange: continuation domain of parameters {[freq1,freq2],[eps1,eps2]},
%           the frequency range should be near the
%           value of natural frequency with index 1 in the mFreq
% outdof:   output for dof in physical domain
% optdof:   amplitude at system level (outdof should be a subset of optdof)
% varargin: [{p0,z0}], ['saveICs'] where {p0,z0} are initial solution
%           guesses and saveICs is a flag saving a point on trajectory as initial
%           condition for numerical integration


m = numel(mFreqs);
assert(numel(modes)==2*m, 'The master subspace is not %dD.',2*m);
nCycle = obj.FRCOptions.nCycle;
%% Checking whether internal resonance indeed happens
if isempty(obj.System.spectrum)
    [~,~,~] = obj.System.linear_spectral_analysis();
end
% eigenvalues Lambda is sorted in descending order of real parts
% positive imaginary part is placed first in each complex pair
lambda = obj.System.spectrum.Lambda(modes);
lambdaRe = real(lambda);
lambdaIm = imag(lambda);
check_spectrum_and_internal_resonance(lambdaRe,lambdaIm,mFreqs);

%% SSM computation of autonomous part
obj.choose_E(modes)
% compute autonomous SSM coefficients
[W_0,R_0] = obj.compute_whisker(order);

% check reduced dynamics (consistency with internal resonance)
[beta,kappa] = check_auto_reduced_dynamics(R_0,order,mFreqs);

% assemble expansion coefficients and indices
wcoeffs = []; cind = []; dind = [];
for k=1:order
    W = W_0(k);
    coeffs = W.coeffs(optdof,:);
    ind = W.ind;
    if ~isempty(coeffs)
        wcoeffs = [wcoeffs coeffs];
        cind = [cind; ind(:,1:2:end-1)];
        dind = [dind; ind(:,2:2:end)];
    end
end

%% SSM computation of non-autonomous part
iNonauto = []; % indices for resonant happens
rNonauto = []; % value of leading order contribution
kNonauto = []; % (pos) kappa indices with resonance
kappa_set= [obj.System.Fext.data.kappa];
kappa_pos = kappa_set(kappa_set>0);
num_kappa = numel(kappa_pos); % number of kappa pairs
for k=1:num_kappa
    kappak = kappa_pos(k);
    idm = find(mFreqs(:)==kappak); % idm could be vector if there are two frequencies are the same
    obj.System.Omega = lambdaIm(2*idm(1)-1);
    [W_1, R_1, Flead] = obj.compute_perturbed_whisker(0,[],[]);
    idk = find(kappa_set==kappak);
    R_10 = R_1(idk).R.coeffs;
    r = R_10(2*idm-1);
    
    iNonauto = [iNonauto; idm];
    rNonauto = [rNonauto; r];  
    kNonauto = [kNonauto; idk];
end

%% create data structure for observables
optData = create_data_for_po_amp(cind,dind,mFreqs,wcoeffs,optdof,obj,Flead);

%% calculation of FRS
if m==1 && strcmp(obj.FRSOptions.method,'level set')
    gammas = compute_gamma(R_0);  
    cal_FRS_via_ana(obj,parRange,gammas,lambda,full(r),outdof,optdof,optData,R_1)
else
    lamd  = struct();
    lamd.lambdaRe = lambdaRe; lamd.lambdaIm = lambdaIm;
    Nonauto = struct();
    Nonauto.iNonauto = iNonauto; Nonauto.rNonauto = rNonauto; Nonauto.kNonauto = kNonauto;
    [fdata,~] = create_reduced_dynamics_data(beta,kappa,lamd,mFreqs,Nonauto,W_0,W_1,order,modes);
    cal_FRS_via_cont_ep(oid,fdata,obj,m,nCycle,scale_state,scale_obs,mFreqs,...
        optData,outdof,optdof,parRange,varargin{:});
end

end