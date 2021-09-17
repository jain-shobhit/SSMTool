function FRC = SSM_cont_ep(obj,type,oid,run,lab,parName,parRange,outdof,varargin)
% SSM_CONT_EP This function performs continuation of equilibrium points of
% slow dynamics. Each equilibrium point corresponds to a periodic orbit in
% the regular time dynamics.
%
% FRC = ODE_CONT_EP(type, OID, RUN, LAB, PARNAME, PARRANGE, OUTDOF, VARARGIN)
%
% type:     BP/ep/SN/HB
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the label of
%           a branch point
% parName:  amp/freq continuation parameter
% parRange: continuation domain of parameter, which should be near the
%           value of natural frequency with index 1 in the mFreq if continuation
%           parameter is freq
% outdof:   dofs for output in physical domain
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: SSM_ISOL2EP SSM_EP2EP

%% continuation of equilibrium points in reduced dynamics
prob = coco_prob();
prob = cocoSet(obj.contOptions,prob);
bif  = false; % true for SN and HB
switch type
    case 'BP'
        prob = ode_BP2ep(prob, '', coco_get_id(run, 'ep'), lab);
    case 'ep'
        prob = ode_ep2ep(prob, '', coco_get_id(run, 'ep'), lab);
    case 'SN'
        prob = ode_ep2SN(prob, '', coco_get_id(run, 'ep'), lab);
        bif  = true;
    case 'HB'
        prob = ode_ep2HB(prob, '', coco_get_id(run, 'ep'), lab);
        bif  = true;
    otherwise
        error('type should be selected from {BP,ep,SN,HB}');
end

% read data 
fdata = coco_get_func_data(prob, 'ep', 'data');
% extract data to fdata.fhan, namely, @(z,p)ode_2mDSSM(z,p,fdata)
odedata = functions(fdata.fhan);
odedata = odedata.workspace{1};
fdata   = odedata.fdata;
m   = numel(fdata.mFreqs);
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
Nonauto = struct();
Nonauto.iNonauto = iNonauto; 
Nonauto.rNonauto = rNonauto; 
Nonauto.kNonauto = kNonauto;


% define monitor functions
[prob, args1, args2] = monitor_states(prob, ispolar, m);

% setup continuation arguments
if ~bif
    switch parName
        case 'freq'
            isomega = true;
            cont_args = [{'om'},args1(:)' ,args2(:)',{'eps'}];
        case 'amp'
            isomega = false;
            cont_args = [{'eps'},args1(:)' ,args2(:)',{'om'}];
        otherwise
            error('Continuation parameter should be freq or amp');
    end
    if strcmp(obj.FRCOptions.sampStyle, 'uniform')
        if isomega
            nOmega = obj.FRCOptions.nPar;
            omSamp = linspace(parRange(1),parRange(2), nOmega);
            prob   = coco_add_event(prob, 'UZ', 'om', omSamp);
        else
            nEpsilon = obj.FRCOptions.nPar;
            epSamp = linspace(parRange(1),parRange(2), nEpsilon);
            prob   = coco_add_event(prob, 'UZ', 'eps', epSamp);
        end        
    end
else
    isomega = true;
    cont_args = [{'om'},{'eps'},args1(:)' ,args2(:)'];  
end

runid = coco_get_id(oid, 'ep');
% print information of continuation run
switch type
    case 'BP'
        fprintf('\n Run=''%s'': Continue equilibria along secondary branch from label %d of run %s.\n', ...
          runid, lab, run);
    case 'ep'
        fprintf('\n Run=''%s'': Continue equilibria from label %d of run %s.\n', ...
          runid, lab, run);
    case 'SN'
        fprintf('\n Run=''%s'': Continue saddle-node equilibria from label %d of run %s.\n', ...
          runid, lab, run);
    case 'HB'
        fprintf('\n Run=''%s'': Continue Hopf equilibria from label %d of run %s.\n', ...
          runid, lab, run);        
end

% coco run
coco(prob, runid, [], 1, cont_args, parRange);


%% extract results of reduced dynamics at sampled frequencies
if ~bif
    FRC = ep_reduced_results(runid,obj.FRCOptions.sampStyle,ispolar,isomega,args1,args2);
else
    FRC = ep_reduced_results(runid,'cocoBD',ispolar,[],args1,args2); 
end

%% FRC in physical domain
FRCdata = struct();        FRCdata.isomega = isomega;
FRCdata.mFreqs  = mFreqs;  FRCdata.order  = order; 
FRCdata.ispolar = ispolar; FRCdata.modes  = modes;
FRC = FRC_reduced_to_full(obj,Nonauto,'ep',FRC,FRCdata,W_0,W_1,outdof,varargin{:});

% Plot Plot FRC in system coordinates
if ~bif
    if isomega
        plot_frc_full(FRC.om,FRC.Znorm_frc,outdof,FRC.Aout_frc,FRC.st,order,'freq','lines',{FRC.SNidx,FRC.HBidx});
    else
        plot_frc_full(FRC.ep,FRC.Znorm_frc,outdof,FRC.Aout_frc,FRC.st,order,'amp','lines',{FRC.SNidx,FRC.HBidx});
    end
else
    if strcmp(type, 'SN')
        plot3_frc_full(FRC.om,FRC.ep,FRC.Znorm_frc,outdof,FRC.Aout_frc,[],order,'lines','g-');
    else
        plot3_frc_full(FRC.om,FRC.ep,FRC.Znorm_frc,outdof,FRC.Aout_frc,[],order,'lines','r-');
    end
end

fdir = fullfile(pwd,'data',runid,'SSMep.mat');
save(fdir, 'FRC');
end