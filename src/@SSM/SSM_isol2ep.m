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
nt = obj.FRCOptions.nt;
nPar = obj.FRCOptions.nPar;
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
flags1 = abs(lambdaRe(1:2:end-1)-lambdaRe(2:2:end))<1e-6*abs(lambdaRe(1:2:end-1)); % same real parts
flags1 = all(flags1);
flags2 = abs(lambdaIm(1:2:end-1)+lambdaIm(2:2:end))<1e-6*abs(lambdaIm(1:2:end-1)); % opposite imag parts
flags2 = all(flags2);
freqs  = lambdaIm(1:2:end-1);
freqso = freqs - dot(freqs,mFreqs(:))*mFreqs(:)/sum(mFreqs.^2);
flags3 = norm(freqso)<0.1*norm(freqs);

assert(flags1, 'Real parts do not follow complex conjugate relation');
assert(flags2, 'Imaginary parts do not follow complex conjugate relation');
assert(flags3, 'Internal resonnace is not detected for given master subspace');

%% SSM computation of autonomous part
obj.choose_E(resonant_modes)
% compute autonomous SSM coefficients
[W_0,R_0] = obj.compute_whisker(order);

% check reduced dynamics to see its consistent with reduced dynamics
beta  = cell(m,1); % coefficients - each cell corresponds to one mode
kappa = cell(m,1); % exponants
Em = eye(m);
for k = 2:order
    R = R_0{k};
    coeffs = R.coeffs;
    ind = R.ind;
    if ~isempty(coeffs)
        for i=1:m
            betai = coeffs(2*i-1,:);
            [~,ki,betai] = find(betai);
            kappai = ind(ki,:);
            % check resonant condition  
            l = kappai(:,1:2:end-1);
            j = kappai(:,2:2:end);
            nk = numel(ki);
            rm = repmat(mFreqs(:)',[nk,1]);
            flagi = dot(l-j-repmat(Em(i,:),[nk,1]),rm,2);
            assert(all(flagi==0), 'Reduced dynamics is not consisent with desired IRs');
            % assemble terms
            beta{i}  = [beta{i} betai];
            kappa{i} = [kappa{i}; kappai];
        end
    end
end

%% SSM computation of non-autonomous part
% 
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
% 
fdata = struct();
fdata.beta  = beta;
fdata.kappa = kappa;
fdata.lamdRe = lambdaRe(1:2:end-1);
fdata.lamdIm = lambdaIm(1:2:end-1);
fdata.mFreqs = mFreqs;
fdata.iNonauto = iNonauto;
fdata.rNonauto = rNonauto;
fdata.kNonauto = kNonauto;
% put W_0 and W_1 in fdata is a bad idea because it will be stored in disk
% for each saved continuation solution. As an alternative, we save W_0 and
% W_1 in disk here under folder data. When needed, they will be loaded into
% memory.
% fdata.W_0   = W_0;
% fdata.W_1   = W_1;
data_dir = fullfile(pwd,'data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
wdir = fullfile(data_dir,'SSM.mat');
SSMcoeffs = struct();
SSMcoeffs.W_0 = W_0;
SSMcoeffs.W_1 = W_1;
save(wdir, 'SSMcoeffs');
fdata.order = order;
fdata.modes = resonant_modes;

ispolar = strcmp(obj.FRCOptions.coordinates, 'polar');
fdata.ispolar = ispolar;
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
% prob = coco_set(prob, 'ode', 'vectorized', false);
prob = cocoSet(obj.contOptions, prob);
p0 = [obj.System.Omega; obj.System.fext.epsilon];
if numel(varargin)>0 && iscell(varargin{1})
    p0 = varargin{1}{1};
    z0 = varargin{1}{2};
    z0 = z0(:);
else
    if ispolar
        z0 = 0.1*ones(2*m,1);
    else
        z0 = zeros(2*m,1);
        % solving linear equations
%         for i=1:numel(iNonauto)
%             id  = iNonauto(i);
%             r   = rNonauto(i);
%             rRe = real(r);
%             rIm = imag(r);
%             ai  = fdata.lamdRe(id);
%             bi  = fdata.lamdIm(id)-mFreqs(id)*p0(1);
%             z0i = [ai -bi;bi ai]\[-rRe;-rIm];
%             z0(2*id-1:2*id) = z0i;
%         end
    end
end
% construct initial guess equilibrium points
switch obj.FRCOptions.initialSolver
    case 'fsolve'
        % fsolve to approximate equilibrium
        fsolveOptions = optimoptions('fsolve','MaxFunctionEvaluations',100000,...
            'MaxIterations',1000000,'FunctionTolerance', 1e-10,...
            'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-10);
        z0 = fsolve(@(z) odefun(z,p0),z0,fsolveOptions);
    case 'forward'
        % forward simulation to approach equilibirum
        tspan = [0 nCycle*2*pi/p0(1)]; %nCycle
        odefw = @(t,z,p) odefun(z,p);
        opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
        [~,y0] = ode45(@(t,y) odefw(t,y,p0), tspan, z0, opts);
        [~,y] = ode45(@(t,y) odefw(t,y,p0), [0 2*pi/p0(1)], y0(end,:));
        [~, warnId] = lastwarn;

        if any(isnan(y(:))) || strcmp(warnId,'MATLAB:ode45:IntegrationTolNotMet')
            warning('Reduced dynamics with IRs in polar form diverges with [0.1 0.1 0.1 0.1]');
        else
            z0 = y(end,:)';
        end
end

if ispolar % regularize initial solution if it is in polar form
    z0(2:2:end) = mod(z0(2:2:end),2*pi); % phase angles in [0,2pi]
    for k=1:m
        if z0(2*k-1)<0
            z0(2*k-1) = -z0(2*k-1);      % positive amplitudes
            z0(2*k) = z0(2*k)+pi;
        end
    end
end

% call ep-toolbox
prob = ode_isol2ep(prob, '', funcs{:}, z0, {'om' 'eps'}, p0);
% define monitor functions to state variables
if ispolar
    rhoargs = cell(m,1);
    thargs  = cell(m,1);
    for k=1:m
        rhoargs{k} = strcat('rho',num2str(k));
        thargs{k}  = strcat('th',num2str(k));
    end
    prob = coco_add_pars(prob, 'radius', 1:2:2*m-1, rhoargs(:)');
    prob = coco_add_pars(prob, 'angle', 2:2:2*m, thargs(:)');
else
    Reargs = cell(m,1);
    Imargs = cell(m,1);
    for k=1:m
        Reargs{k} = strcat('Rez',num2str(k));
        Imargs{k} = strcat('Imz',num2str(k));
    end
    prob = coco_add_pars(prob, 'realParts', 1:2:2*m-1, Reargs(:)');
    prob = coco_add_pars(prob, 'imagParts', 2:2:2*m, Imargs(:)');
end    

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
    if ispolar
        cont_args = [{'om'},rhoargs(:)' ,thargs(:)',{'eps'}];
    else
        cont_args = [{'om'},Reargs(:)' ,Imargs(:)',{'eps'}];
    end
else
    if ispolar
        cont_args = [{'eps'},rhoargs(:)' ,thargs(:)',{'om'}];
    else
        cont_args = [{'eps'},Reargs(:)' ,Imargs(:)',{'om'}];
    end
end
    

coco(prob, runid, [], 1, cont_args, parRange);

%% extract results of reduced dynamics at sampled frequencies
if ispolar
    FRC = ep_reduced_results(runid,obj.FRCOptions.sampStyle,ispolar,isomega,rhoargs,thargs);
else
    FRC = ep_reduced_results(runid,obj.FRCOptions.sampStyle,ispolar,isomega,Reargs,Imargs);    
end

%% FRC in physical domain
Zout_frc = [];
Znorm_frc = [];
Aout_frc = [];

% flag for saving ICs (used for numerical integration)
if numel(varargin)==2
    saveICflag = strcmp(varargin{2},'saveICs');
elseif numel(varargin)==1
    saveICflag = strcmp(varargin{1},'saveICs');
else
    saveICflag = false;
end
if saveICflag
    Z0_frc = []; % initial state
end
timeFRCPhysicsDomain = tic;

% Loop around a resonant mode
om   = FRC.om;
epsf = FRC.ep;
for j = 1:numel(om)
    % compute non-autonomous SSM coefficients
    obj.System.Omega = om(j);
    if isomega
        if obj.Options.contribNonAuto
            [W_1, R_1] = obj.compute_perturbed_whisker(order);

            R_10 = R_1{1}.coeffs;
            assert(~isempty(R_10), 'Near resonance does not occur, you may tune tol');
            f = R_10((kNonauto-1)*2*m+2*iNonauto-1);

            assert(norm(f-rNonauto)<1e-3*norm(f), 'inner resonance assumption does not hold');

            fprintf('the forcing frequency %.4d is nearly resonant with the eigenvalue %.4d + i%.4d\n', om(j), real(lambda(1)),imag(lambda(1)))
        else
            W_1 = [];
        end
    end
    % Forced response in Physical Coordinates
    state = FRC.z(j,:);
    [Aout, Zout, z_norm, Zic] = compute_full_response_2mD_ReIm(W_0, W_1, state, epsf(j), nt, mFreqs, outdof);

    % collect output in array
    Aout_frc = [Aout_frc; Aout];
    Zout_frc = [Zout_frc; Zout];
    Znorm_frc = [Znorm_frc; z_norm];

    if saveICflag
        Z0_frc = [Z0_frc; Zic]; % initial state
    end 
end
%% 
% Record output
FRC.Aout_frc  = Aout_frc;
FRC.Zout_frc  = Zout_frc;
FRC.Znorm_frc = Znorm_frc;
if saveICflag
    FRC.Z0_frc = Z0_frc; % initial state
end 
FRC.timeFRCPhysicsDomain = toc(timeFRCPhysicsDomain);
FRC.SSMorder   = order;
FRC.SSMnonAuto = obj.Options.contribNonAuto;
FRC.SSMispolar = ispolar;

% Plot Plot FRC in system coordinates
if isomega
    plot_frc_full(om,Znorm_frc,outdof,Aout_frc,FRC.st,order,'freq','lines',{FRC.SNidx,FRC.HBidx});
else
    plot_frc_full(epsf,Znorm_frc,outdof,Aout_frc,FRC.st,order,'amp','lines',{FRC.SNidx,FRC.HBidx});
end

varargout{1} = FRC;
fdir = fullfile(data_dir,runid,'SSMep.mat');
save(fdir, 'FRC');
end
