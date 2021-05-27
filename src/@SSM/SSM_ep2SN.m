function varargout = SSM_ep2SN(obj,oid,run,lab,parRange,outdof,varargin)
% SSM_EP2SN This function performs continuation of saddle-node (SN)
% equilibirium points of slow dynamics. SN bifurcation is of codimension
% one and hence two parameters are free to vary to yield an one-dimensional
% manifold of SN points. Each SN point corresponds to a SN bifurcation
% periodic orbit in the regular time dynamics. The continuation here starts
% from a saved solution, which is a SN point.
%
% FRCIRS = SSM_EP2SN(OBJ,OID,RUN,LAB,PARRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of current continuation
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the
%           label of a saddle-node point
% parRange: continuation domain of parameters. It is of the form
%           {[om1,om2],[f1,f2]}, where [om1,om2] and [f1,f2] specify the
%           continuation domain of excitation frequency and amplitude
%           respectively. You can give empty array and then no domain is
%           specified, e.g., {[],[f1,f2]} only presents the domain of
%           forcing amplitude
% outdof:   output for dof in physical domain
%
% See also: SSM_ISOL2EP, SSM_EP2EP, SSM_BP2EP

%% continuation of equilibrium points in reduced dynamics
prob = coco_prob();
prob = cocoSet(obj.contOptions,prob);
prob = ode_ep2SN(prob, '', coco_get_id(run, 'ep'), lab);

% read data 
fdata = coco_get_func_data(prob, 'ep', 'data');
% extract data to fdata.fhan, namely, @(z,p)ode_2mDSSM(z,p,fdata)
odedata = functions(fdata.fhan);
odedata = odedata.workspace{1};
fdata   = odedata.fdata;
m   = numel(fdata.mFreqs);
% W_0 = fdata.W_0;
% W_1 = fdata.W_1;
order    = fdata.order;
ispolar  = fdata.ispolar;
iNonauto = fdata.iNonauto;
rNonauto = fdata.rNonauto;
kNonauto = fdata.kNonauto;
modes    = fdata.modes;
mFreqs   = fdata.mFreqs;
wdir     = fullfile(pwd,'data','SSM.mat');
SSMcoeffs = load(wdir);
SSMcoeffs = SSMcoeffs.SSMcoeffs;
W_0 = SSMcoeffs.W_0;
W_1 = SSMcoeffs.W_1;
clear('SSMcoeffs');

% define monitor functions
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

% setup continuation arguments
runid = coco_get_id(oid, 'ep');
fprintf('\n Run=''%s'': Continue saddle-node equilibria.\n', ...
  runid);
if ispolar
    cont_args = [{'om'},{'eps'},rhoargs(:)',thargs(:)'];
else
    cont_args = [{'om'},{'eps'},Reargs(:)',Imargs(:)'];
end
assert(iscell(parRange),'the parRange should be a cell with two entries');

% coco run
coco(prob, runid, [], 1, cont_args, parRange);


%% extract results of reduced dynamics at sampled frequencies
if ispolar
    FRC = ep_reduced_results(runid,'cocoBD',ispolar,[],rhoargs,thargs);
else
    FRC = ep_reduced_results(runid,'cocoBD',ispolar,[],Reargs,Imargs);    
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
nt   = obj.FRCOptions.nt;
for j = 1:numel(om)
    % compute non-autonomous SSM coefficients
    obj.System.Omega = om(j);
    if obj.Options.contribNonAuto
        if isempty(obj.E)
            obj.choose_E(modes);
        end

        [W_1, R_1] = obj.compute_perturbed_whisker(order);

        R_10 = R_1{1}.coeffs;
        assert(~isempty(R_10), 'Near resonance does not occur, you may tune tol');
        f = R_10((kNonauto-1)*2*m+2*iNonauto-1);

        assert(norm(f-rNonauto)<1e-3*norm(f), 'inner resonance assumption does not hold');

        fprintf('the forcing frequency %.4d is nearly resonant with the eigenvalue %.4d + i%.4d\n',...
            om(j), fdata.lamdRe(1),fdata.lamdIm(1))
    else
        W_1 = [];
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
plot3_frc_full(om,epsf,Znorm_frc,outdof,Aout_frc,[],order,'lines','g-');

varargout{1} = FRC;
fdir = fullfile(pwd,'data',runid,'SSMep.mat');
save(fdir, 'FRC');
end