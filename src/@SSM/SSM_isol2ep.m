function varargout = SSM_isol2ep(obj,oid,resonant_modes,order,mFreqs,parName,parRange,outdof,varargin)
% SSM_ISOL2EP This function performs continuation of equilibrium points of
% slow dynamics. Each equilibirum point corresponds to a periodic orbit in
% the regular time dynamics. The continuation here starts from the guess of
% initial solution.
%
% FRC = SSM_ISOL2EP(OBJ,OID,RESONANT_MODES,ORDER,MFREQS,PARNAME,PARRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of continuation
% resonant_modes:    master subspace
% order:    expansion order of SSM
% mFreqs:   internal resonance relation vector
% parName:  amp/freq continuation parameter
% parRange: continuation domain of parameter, which should be near the
%           value of natural frequency with index 1 in the mFreq if the continuation
%           parameter is freq
% outdof:   output for dof in physical domain
% varargin: [{p0,z0}], ['saveICs'] where {p0,z0} are initial solution
%           guesses and saveICs is a flag saving a point on trajectory as initial
%           condition for numerical integration


m = numel(mFreqs);
assert(numel(resonant_modes)==2*m, 'The master subspace is not %dD.',2*m);
nPar   = obj.FRCOptions.nPar;
nCycle = obj.FRCOptions.nCycle;
%% Checking whether internal resonance indeed happens
if isempty(obj.System.spectrum)
    [~,~,~] = obj.System.linear_spectral_analysis();
end
% eigenvalues Lambda is sorted in descending order of real parts
% positive imaginary part is placed first in each complex pair

lambda = obj.System.spectrum.Lambda(resonant_modes);
lambdaRe = real(lambda);
lambdaIm = imag(lambda);
check_spectrum_and_internal_resonance(lambdaRe,lambdaIm,mFreqs);

%% SSM computation of autonomous part
obj.choose_E(resonant_modes)
% compute autonomous SSM coefficients
[W_0,R_0] = obj.compute_whisker(order);

% check reduced dynamics (consistent with expected internal resonance)
[beta,kappa] = check_auto_reduced_dynamics(R_0,order,mFreqs);

%% SSM computation of non-autonomous part
iNonauto = []; % indices for resonant happens
rNonauto = []; % value of leading order contribution
kNonauto = []; % (pos) kappa indices with resonance
kappa_set= obj.System.Fext.kappas; % each row corresponds to one kappa
kappa_pos = kappa_set(kappa_set>0);
num_kappa = numel(kappa_pos); % number of kappa pairs
for k=1:num_kappa
    kappak = kappa_pos(k);
    idm = find(mFreqs(:)==kappak); % idm could be vector if there are two frequencies are the same
    obj.System.Omega = lambdaIm(2*idm(1)-1);
    
    [W_1, R_1] = obj.compute_perturbed_whisker(order);

    R_10 = R_1{1}.coeffs;
    idk = find(kappa_set==kappak);
    r = R_10(2*idm-1,idk);
    
    iNonauto = [iNonauto; idm];
    rNonauto = [rNonauto; r];  
    kNonauto = [kNonauto; idk];
end
%% Construct COCO-compatible vector field 
% create data to vector field 
lamd  = struct(); 
lamd.lambdaRe = lambdaRe; lamd.lambdaIm = lambdaIm;
Nonauto = struct();
Nonauto.iNonauto = iNonauto; Nonauto.rNonauto = rNonauto; Nonauto.kNonauto = kNonauto;
[fdata,data_dir] = create_reduced_dynamics_data(beta,kappa,lamd,mFreqs,Nonauto,W_0,W_1,order,resonant_modes);

ispolar = strcmp(obj.FRCOptions.coordinates, 'polar');
fdata.ispolar = ispolar;
fdata.isbaseForce = obj.System.Options.BaseExcitation;
if ispolar
    odefun = @(z,p) ode_2mDSSM_polar(z,p,fdata);
else
    odefun = @(z,p) ode_2mDSSM_cartesian(z,p,fdata);
end
funcs  = {odefun};
% 
%% continuation of reduced dynamics w.r.t. parName
%
prob = coco_prob();
prob = cocoSet(obj.contOptions, prob);
% construct initial solution
if obj.System.order==2
    p0 = [obj.System.Omega; obj.System.fext.epsilon];
else
    p0 = [obj.System.Omega; obj.System.Fext.epsilon];
end
[p0,z0] = initial_fixed_point(p0,obj.FRCOptions.initialSolver,ispolar,...
    odefun,nCycle,m,varargin{:});

% call ep-toolbox
prob = ode_isol2ep(prob, '', funcs{:}, z0, {'om' 'eps'}, p0);
% define monitor functions to state variables
[prob, args1, args2] = monitor_states(prob, ispolar, m);

switch parName
    case 'freq'
        isomega = true;
    case 'amp'
        isomega = false;
    otherwise
        error('Continuation parameter should be freq or amp');
end
if strcmp(obj.FRCOptions.sampStyle, 'uniform')
    if isomega
        omSamp = linspace(parRange(1),parRange(2), nPar);
        prob   = coco_add_event(prob, 'UZ', 'om', omSamp);
    else
        epSamp = linspace(parRange(1),parRange(2), nPar);
        prob   = coco_add_event(prob, 'UZ', 'eps', epSamp);
    end        
end

runid = coco_get_id(oid, 'ep');

fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  runid);

if isomega
    cont_args = [{'om'},args1(:)' ,args2(:)',{'eps'}];
else
    cont_args = [{'eps'},args1(:)' ,args2(:)',{'om'}];
end
    
coco(prob, runid, [], 1, cont_args, parRange);

%% extract results of reduced dynamics at sampled frequencies
FRC = ep_reduced_results(runid,obj.FRCOptions.sampStyle,ispolar,isomega,args1,args2);

%% FRC in physical domain
FRCdata = struct();        FRCdata.isomega = isomega; 
FRCdata.mFreqs  = mFreqs;  FRCdata.order  = order; 
FRCdata.ispolar = ispolar; FRCdata.modes  = resonant_modes;
FRC = FRC_reduced_to_full(obj,Nonauto,'ep',FRC,FRCdata,W_0,W_1,outdof,varargin{:});

% Plot Plot FRC in system coordinates
if isomega
    plot_frc_full(FRC.om,FRC.Znorm_frc,outdof,FRC.Aout_frc,FRC.st,order,'freq','lines',{FRC.SNidx,FRC.HBidx});
else
    plot_frc_full(FRC.ep,FRC.Znorm_frc,outdof,FRC.Aout_frc,FRC.st,order,'amp','lines',{FRC.SNidx,FRC.HBidx});
end

varargout{1} = FRC;
fdir = fullfile(data_dir,runid,'SSMep.mat');
save(fdir, 'FRC');
end
