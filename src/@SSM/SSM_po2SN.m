function SSM_po2SN(obj,oid,run,lab,parRange,outdof,varargin)
% SSM_PO2SN This function performs continuation of saddle-node (SN)
% periodic orbits of slow dynamics. Each periodic orbit corresponds to a
% torus (quasi-periodic) response in regular time dynamics. The
% continuation here starts from a saved SN solution.
%
% FRCIRS = SSM_PO2SN(OBJ,OID,RUN,LAB,PARRANGE,OUTDOF,VARARGIN)
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
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: SSM_HB2PO, SSM_PO2PO, SSM_BP2PO

%% continuation of equilibrium points in reduced dynamics
prob = coco_prob();
prob = cocoSet(obj.contOptions,prob);
prob = ode_po2SN(prob, '', coco_get_id(run, 'po'), lab);

% read data 
fdata = coco_get_func_data(prob, 'po.orb.coll', 'data');
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
dim      = 2*m;
wdir     = fullfile(pwd,'data','SSM.mat');
SSMcoeffs = load(wdir);
SSMcoeffs = SSMcoeffs.SSMcoeffs;
W_0 = SSMcoeffs.W_0;
W_1 = SSMcoeffs.W_1;
clear('SSMcoeffs');

% setup continuation arguments
cont_args = {'om', 'eps', 'po.period'};
runid = coco_get_id(oid, 'po');
fprintf('\n Run=''%s'': Continue SN periodic orbits from label %d of run %s.\n', ...
  runid, lab, run);

% coco run
coco(prob, runid, [], 1, cont_args, parRange);

%% extract results of reduced dynamics at sampled frequencies
FRC = po_reduced_results(runid,ispolar,[],mFreqs,obj.FRCOptions.nt,false);

%% map torus back to physical system
fprintf('\n FRCs from =''%s'': generating torus in physical domain.\n', ...
    runid);
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
nlab = numel(FRC.lab);
zTr  = cell(nlab,1);
om   = FRC.om;
epsf = FRC.ep;
nSeg = FRC.nSeg;
noutdof = numel(outdof);
for j = 1:nlab
    % compute non-autonomous SSM coefficients
    obj.System.Omega = om(j);
    switch obj.Options.contribNonAuto
        case 'keep'
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
        case 'delete'
            W_1 = [];
    end
    % Forced response in physical Coordinates
    qTrj = FRC.qTr{j};
    tTrj = FRC.tTr{j};
    nt   = numel(tTrj);
    numSegs = nSeg(j);
    Zout_frc = zeros(nt,noutdof,numSegs);
    for k=1:numSegs
        xbp = qTrj(:,:,k);
        xbp = xbp';
        x_comp = xbp(1:2:end-1,:)+1i*xbp(2:2:end,:);
        state  = zeros(dim, nt); % state
        state(1:2:end-1,:) = x_comp;
        state(2:2:end,:)   = conj(x_comp);
        
        [~, Zout, ~,Zic] = compute_full_response_traj(W_0, W_1, epsf(j), tTrj, state, om(j), outdof);
        Zout_frc(:,:,k) = Zout';
    end
    
    zTr{j} = Zout_frc;

    if saveICflag
        Z0_frc = [Z0_frc Zic]; % initial state
    end  
end
%% 
% Record output
FRC.zTr  = zTr;
if saveICflag
    FRC.Z0_frc = Z0_frc; % initial state
end
FRC.timeFRCPhysicsDomain = toc(timeFRCPhysicsDomain);
FRC.SSMorder   = order;
FRC.SSMnonAuto = obj.Options.contribNonAuto;
FRC.SSMispolar = ispolar;

fdir = fullfile(pwd,'data',runid,'SSMpo.mat');
save(fdir, 'FRC');
end
