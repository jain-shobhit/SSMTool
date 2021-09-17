function varargout = SSM_cont_tor(obj,type,oid,run,lab,parName,parRange,outdof,varargin)
% SSM_CONT_TOR This function performs continuation of 2D tori of slow
% dynamics. Each 2D tori corresponds to a 3D torus (quasi-periodic)
% response in regular time dynamics. 
%
% FRCIRS = SSM_TR2TOR(OBJ,OID,RUN,LAB,PARNAME,PARRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of current continuation
% type:     TR/tor/BP
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the
%           label of a TR/BP point if type=TR/BP
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
% See also: SSM_PO2TR, SSM_TOR2PTOR, SSM_CONT_TOR

%% continuation of equilibrium points in reduced dynamics
prob = coco_prob();
prob = coco_set(prob, 'tor', 'autonomous', true, 'nOmega', 0);
prob = cocoSet(obj.contOptions,prob);
runid = coco_get_id(oid, 'tor');
switch type
    case 'TR'
        prob = ode_TR2tor(prob, '', coco_get_id(run, 'po'), lab,...
            obj.FRCOptions.torNumSegs, obj.FRCOptions.torRotDiret,...
            obj.FRCOptions.torPurtb);
        fprintf('\n Run=''%s'': Continue tori born from TR point with label %d of run %s.\n', ...
          runid, lab, run);
    case 'tor'
        prob = ode_tor2tor(prob, '', coco_get_id(run, 'tor'), lab);
        fprintf('\n Run=''%s'': Continue tori from label %d of run %s.\n', ...
          runid, lab, run);        
    case 'BP'
        prob = ode_BP2tor(prob, '', coco_get_id(run, 'tor'), lab);
        fprintf('\n Run=''%s'': Continue tori along secondary branch from label %d of run %s.\n', ...
          runid, lab, run);        
    otherwise
        error('type should be selected from {isol,HB,BP,PO,SN,HB}');
end        

% read data 
fdata = coco_get_func_data(prob, 'tor.bvp.seg1.coll', 'data');
% extract data to fdata.fhan, namely, @(z,p)ode_2mDSSM(z,p,fdata)
odedata = functions(fdata.fhan);
odedata = odedata.workspace{1};
odedata = functions(odedata.args.fhan);
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

% coco run
coco(prob, runid, [], 1, cont_args, parRange);

%% extract results of reduced dynamics at sampled frequencies
FRC = tor_reduced_results(runid,ispolar,isomega,mFreqs,obj.FRCOptions.nt);

%% map torus back to physical system
fprintf('\n FRCs from =''%s'': generating torus in physical domain.\n', ...
    runid);
FRCdata = struct();        FRCdata.isomega = isomega;
FRCdata.mFreqs  = mFreqs;  FRCdata.order  = order; 
FRCdata.ispolar = ispolar; FRCdata.modes  = modes;   FRCdata.dim = dim;
FRC = FRC_reduced_to_full(obj,Nonauto,'tor',FRC,FRCdata,W_0,W_1,outdof,varargin{:});
varargout{1} = FRC;
% write results to hard driver
fdir = fullfile(pwd,'data',runid,'SSMtor.mat');
save(fdir, 'FRC');
end
