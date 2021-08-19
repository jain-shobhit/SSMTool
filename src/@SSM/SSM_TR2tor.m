function SSM_TR2tor(obj,oid,run,lab,parName,parRange,outdof,varargin)
% SSM_TR2tor This function performs continuation of 2D tori of slow
% dynamics. Each 2D tori corresponds to a 3D torus (quasi-periodic)
% response in regular time dynamics. The continuation here starts by
% switching to branch of tori at torus (Neimark-Sacker) bifurcation point
% found in the continuation of periodic orbits. 
%
% FRCIRS = SSM_TR2TOR(OBJ,OID,RUN,LAB,PARNAME,PARRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of current continuation
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the
%           label of a torus point
% parName:  freq-amp/amp/freq continuation parameter. In the case of
%           freq-amp, both two parameters are free to change and hence
%           varrho is fixed. In contrast, varrho is fixed if there is only
%           one continuation parameter
% parRange: continuation domain of parameter, which is given in the form of
%           {[om1,om2],[f1,f2]} for freq-amp, [f1,f2] for amp and [om1,om2]
%           for freq. It is noted that the domain [om1,om2] should contain
%           the natural frequency with index 1 in the mFreq vector
% outdof:   output for dof in physical domain
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: SSM_PO2TR, SSM_TOR2PTOR

%% continuation of equilibrium points in reduced dynamics
prob = coco_prob();
prob = coco_set(prob, 'tor', 'autonomous', true, 'nOmega', 0);
prob = cocoSet(obj.contOptions,prob);
prob = ode_TR2tor(prob, '', coco_get_id(run, 'po'), lab,...
    obj.FRCOptions.torNumSegs, obj.FRCOptions.torRotDiret,...
    obj.FRCOptions.torPurtb);

% read data 
fdata = coco_get_func_data(prob, 'tor.bvp.seg1.coll', 'data');
% extract data to fdata.fhan, namely, @(z,p)ode_2mDSSM(z,p,fdata)
odedata = functions(fdata.fhan);
odedata = odedata.workspace{1};
odedata = functions(odedata.args.fhan);
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
switch parName
    case 'freq-amp'
        % varrho is fixed
        cont_args = {'om', 'eps', 'om1' ,'om2', 'varrho'};
    case 'freq'
        % varrho is free to vary
        isomega = true;
        cont_args = {'om', 'varrho', 'om1', 'om2', 'eps'};
    case 'amp'
        % varrho is free to vary
        isomega = false;
        cont_args = {'eps', 'varrho', 'om1', 'om2', 'om'};
    otherwise
        error('Continuation parameter should be freq or amp');
end

runid = coco_get_id(oid, 'tor');
fprintf('\n Run=''%s'': Continue tori born from TR point.\n', ...
  runid);

% coco run
coco(prob, runid, [], 1, cont_args, parRange);


%% extract results of reduced dynamics at sampled frequencies
FRC = tor_reduced_results(runid,ispolar,isomega,mFreqs,obj.FRCOptions.nt);

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
    if isomega
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

fdir = fullfile(pwd,'data',runid,'SSMtor.mat');
save(fdir, 'FRC');
end
