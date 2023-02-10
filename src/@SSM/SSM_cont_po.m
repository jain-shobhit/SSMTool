function varargout = SSM_cont_po(obj,type,oid,run,lab,initsol,parName,parRange,outdof,varargin)
% SSM_CONT_PO This function performs continuation of periodic orbits of slow
% dynamics. Each periodic orbit corresponds to a torus (quasi-periodic)
% response in regular time dynamics.
%
% FRCIRS = SSM_CONT_PO(OBJ,TYPE,OID,RUN,LAB,INITSOL,PARNAME,PARRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of current continuation
% type:     isol/HB/BP/po/SN/PD/TR
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the
%           label of a HB/BP/SN/PD/TR point if type=HB/BP/SN/PD/TR
% initsol:  initial solution (type=isol), which is a structure (t,x,eps,om)
% parName:  amp/freq continuation parameter
% parRange: continuation domain of parameter, which should be near the
%           value of natural frequency with index 1 in the mFreq if continuation
%           parameter is freq
% outdof:   output for dof in physical domain
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: SSM_EP2HB, SSM_PO2PO

%% continuation of equilibrium points in reduced dynamics
prob  = coco_prob();
prob  = cocoSet(obj.contOptions,prob);
bif   = false; % true for SN/PD/TR
runid = coco_get_id(oid, 'po');
switch type
    case 'isol'
        % continuation from an initial periodic solution guess
        [~, data] = ep_read_solution('', coco_get_id(run, 'ep'), lab);
        t0 = initsol.t;
        x0 = initsol.x;
        p0 = [initsol.om;initsol.eps];
        prob = ode_isol2po(prob, '', data.fhan, data.dfdxhan, data.dfdphan, ...
          t0, x0, data.pnames, p0);
        fprintf('\n Run=''%s'': Continue periodic orbits with initial solution.\n', ...
          runid);
    case 'HB'
        % continuation from a Hopf bifurcation fixed point
        prob = ode_HB2po(prob, '', coco_get_id(run, 'ep'), lab);   
        fprintf('\n Run=''%s'': Continue periodic orbits born from a HB point with label %d of run %s.\n', ...
          runid, lab, run);
    case 'BP'
        % continuation along a secondary branch passing a branch point
        prob = ode_BP2po(prob, '', coco_get_id(run, 'po'), lab);  
        fprintf('\n Run=''%s'': Continue periodic orbits along secondary branch of solution with label %d of run %s.\n', ...
          runid, lab, run);        
    case 'po'
        % continuation from a saved periodic solution
        prob = ode_po2po(prob, '', coco_get_id(run, 'po'), lab);  
        fprintf('\n Run=''%s'': Continue periodic orbits born from saved solution with label %d of run %s.\n', ...
          runid, lab, run);        
    case 'SN'
        % continuation of saddle-node bifurcation periodic orbit
        prob = ode_po2SN(prob, '', coco_get_id(run, 'po'), lab);
        bif  = true;
        fprintf('\n Run=''%s'': Continue SN periodic orbits from label %d of run %s.\n', ...
          runid, lab, run);    
    case 'PD'
        % continuation of period-doubling bifurcation periodic orbit
        prob = ode_po2PD(prob, '', coco_get_id(run, 'po'), lab);
        bif  = true;
        fprintf('\n Run=''%s'': Continue PD periodic orbits from label %d of run %s.\n', ...
          runid, lab, run);        
    case 'TR'
        % continuation of torus bifurcation periodic orbit
        prob = ode_po2TR(prob, '', coco_get_id(run, 'po'), lab);
        bif  = true;
        fprintf('\n Run=''%s'': Continue TR periodic orbits from label %d of run %s.\n', ...
          runid, lab, run); 
    case 'Tinf'
        prob = ode_po2po(prob, '', coco_get_id(run, 'po'), lab);
        Tidx = coco_get_func_data(prob,'po.period','uidx');
        prob = coco_add_pars(prob,'fixT',Tidx,'T');
        bif = true;
    otherwise
        error('type should be selected from {isol,HB,BP,PO,SN,HB}');
end

% read data 
fdata = coco_get_func_data(prob, 'po.orb.coll', 'data');
% extract data to fdata.fhan, namely, @(z,p)ode_2mDSSM(z,p,fdata)
odedata = functions(fdata.fhan);
odedata = odedata.workspace{1};
fdata   = odedata.fdata;
m       = numel(fdata.mFreqs);
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
Nonauto = struct();
Nonauto.iNonauto = iNonauto; 
Nonauto.rNonauto = rNonauto; 
Nonauto.kNonauto = kNonauto;

% setup continuation arguments
isuniform = false;
if ~bif
    switch parName
        case 'freq'
            isomega = true;
            cont_args = {'om', 'po.period', 'eps'};
        case 'amp'
            isomega = false;
            cont_args = {'eps', 'po.period', 'om'};
        otherwise
            error('Continuation parameter should be freq or amp');
    end

    if ~isempty(obj.FRCOptions.parSamps)
        if isomega
            prob = coco_add_event(prob, 'PS','om',obj.FRCOptions.parSamps);
        else
            prob = coco_add_event(prob, 'PS','eps',obj.FRCOptions.parSamps);
        end
    end

    if strcmp(obj.FRCOptions.sampStyle, 'uniform')
        isuniform = true;
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
    cont_args = {'om', 'eps', 'po.period'};
    if strcmp(type,'Tinf')
        isomega = true;
    else
        isomega = [];
    end
end

% coco run
coco(prob, runid, [], 1, cont_args, parRange);

%% extract results of reduced dynamics at sampled frequencies
FRC = po_reduced_results(runid,ispolar,isomega,mFreqs,obj.FRCOptions.nt,isuniform);

%% map torus back to physical system
if isempty(isomega) 
    isomega = true; 
end
fprintf('\n FRCs from =''%s'': generating torus in physical domain.\n', ...
    runid);
FRCdata = struct();        FRCdata.isomega = isomega;
FRCdata.mFreqs  = mFreqs;  FRCdata.order  = order; 
FRCdata.ispolar = ispolar; FRCdata.modes  = modes;   FRCdata.dim = dim;
FRC = FRC_reduced_to_full(obj,Nonauto,'po',FRC,FRCdata,W_0,W_1,outdof,varargin{:});
varargout{1} = FRC;
% write results to hard driver
fdir = fullfile(pwd,'data',runid,'SSMpo.mat');
save(fdir, 'FRC');
end
