function varargout = FRC_cont_po(obj,oid,resModes,order,parRange)
% FRC_CONT_PO This function performs continuation of periodic orbits of
% the reduced dynamics with respect to the forcing frequency for 2 dimensional SSMs. 
%
% FRC = FRC_CONT_PO(OBJ,OID,RESMODES,ORDER,PARRANGE)
%
% oid:          runid of continuation
% resModes:     master subspace
% order:        expansion order of SSM
% parRange:     continuation domain of parameter, which should be near the
%               value of natural frequency or around 2:1 resonance of it
%
% varargout:    cell array that contains FRC struct
%
% See also: FRC_CONT_EP, FRC_LEVEL_SET, EXTRACT_FRC

% get options
[nt,nCycle,sampStyle,outdof,saveIC,coordinates,periodsRatio,p0,z0] = ...
    deal(obj.FRCOptions.nt,  obj.FRCOptions.nCycle, obj.FRCOptions.sampStyle,  ...
    obj.FRCOptions.outdof, obj.FRCOptions.saveIC, obj.FRCOptions.coordinates, ...
    obj.FRCOptions.periodsRatio, obj.FRCOptions.p0, obj.FRCOptions.z0);

% Check if method is applicable
dimModes = numel(resModes);
assert(dimModes==2,'continuation_po method for FRC extraction is implemented only for 2-dimensional SSMs. Please use the continuation method')

%% SSM computation of autonomous part
obj.choose_E(resModes)
[W,R] = obj.compute_whisker(order);


%% SSM computation of nonautonomous part

if isempty(p0)
    obj.System.Omega = parRange(1);
else
    obj.System.Omega = p0(1);
end

%% Construct COCO-compatible vector field and Jacobians

ispolar = strcmp(coordinates, 'polar');

if ispolar
    error('ODE in polar coordinates not implemented for method continuation po, please change to cartesian')
else
    
    if obj.FRCOptions.omDepNonAuto % Assume strong dependence of coefficients on Omega
        
        fdata   = struct('order', order,'R',R,'W',W);
        odefun    = @(t,x,p) ode_2DSSM_cartesian(t,x,p,fdata,obj);
        odefun_dx = @(t,x,p) ode_2DSSM_cartesian_DFDX(t,x,p,fdata,obj);
        odefun_dp = @(t,x,p) ode_2DSSM_cartesian_DFDP(t,x,p,fdata,obj);
    else
                
        % Compute Reduced dynamics and sensitivity coefficients
        obj.System.Omega = obj.FRCOptions.omDepNonAutoVal;
        [X, S] = obj.compute_perturbed_whisker(order-1,W,R);
        % [~,DS] = obj.compute_sensitivity_coefficients(order-1,W,R,X);

        
        fdata   = struct('order', order,'R',R,'S',S); %,'DS',DS);
        odefun    = @(t,x,p) ode_2DSSM_cartesian_fixROM(t,x,p,fdata);
        odefun_dx = @(t,x,p) ode_2DSSM_cartesian_fixROM_DFDX(t,x,p,fdata);
        odefun_dp = @(t,x,p) ode_2DSSM_cartesian_fixROM_DFDP(t,x,p,fdata,obj.Options.contribNonAuto);
    end
    
    funcs  = {odefun,odefun_dx,odefun_dp};
end

%% initial solution by forward simulation
% ode45 is used here. Integration option may be added in future
if isempty(p0)
    p0 = [parRange(1); obj.System.Fext.epsilon];
end

if isempty(z0)
    if ispolar
        z0 = 0.1*ones(2,1);
    else
        z0 = zeros(2,1);
    end
end

% construct initial guess periodic orbit
[z0,t0] = get_initial_sol(z0,p0,odefun,nCycle,obj.System.Omega,ispolar,periodsRatio);

%% Build continuation problem
% setup coco
prob = coco_prob();
prob = cocoSet(obj.contOptions, prob); %set default options 
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', false);


coll_args = {funcs{:}, t0, z0, {'om','eps'}, p0};   %#ok<CCAT>
prob = ode_isol2po(prob, '', coll_args{:});

% store floquet multipliers in bd
prob = po_mult_add(prob,'po.orb'); 

% Fix response period
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
glueperiod    = @(prob, data, u)  deal(data,u(1,:) - periodsRatio* (2*pi / u(2,:)));
glueperiod_du = @(prob, data, u)  deal(data,[1,  periodsRatio*(2*pi) / (u(2)^2)]);

prob = coco_add_func(prob, 'OmegaT', glueperiod,glueperiod_du,  [],'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);


%% Monitor function for full system response
% Convert reduced results to full while continuing to avoid overhead of SSM
% computations
ampData = struct('W',W,'nt',nt,'outdof',outdof,'saveIC',saveIC);
numoutdof   = numel(outdof);

% Set up problem to either save IC of every p.o. or not.
ampNames    = cell(1, numoutdof+1);
ampNames{end} = 'Znorm';
for k = 1:numoutdof
    ampNames{k} = strcat('amp',num2str(outdof(k)));
end


ampfunc = @(prob,data,u) full_amplitude(prob,data,u,obj);
prob = coco_add_func(prob, 'amp', ampfunc, ampData, 'regular', ampNames,...
    'uidx', uidx([maps.x0_idx, maps.p_idx]));

% monitor function for saving the initial condition
if saveIC
    prob = coco_add_slot(prob, 'x0tobd', @init_cond, data, 'bddat');
end
%% Start continuation
cont_args = {[{'om'},{'po.period'},ampNames(:)'],parRange};

% Release parameters for continuation
runid = coco_get_id(oid, 'po');
fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
    runid);


bd = coco(prob, runid, [], cont_args{:});

%% FRC in physical domain
% extract data from bifurcation data of primary branch
FRC{1} = po_full_results(runid,sampStyle,ispolar,ampNames,saveIC,'primary',[],'plot-off'); 
%% Continue secondary branches if there are any
if  obj.FRCOptions.branchSwitch 
    
    labs = coco_bd_labs(bd, 'BP');
    
    ampCell = {ampData,ampNames,ampfunc};
    periodCell = {glueperiod,glueperiod_du};
    ii = 2;
    for i = 1:numel(labs)
        lab = labs(i);
        runid_i = strcat(runid,'_BP_',num2str(i));
        
        % Continue branchpoint
        continuation_BP(runid,runid_i,lab,obj.contOptions,cont_args,periodCell,ampCell, saveIC);
        
        % Add data of secondary branches to the FRC struct
        FRC{ii} = po_full_results(runid_i,sampStyle,ispolar,ampNames,saveIC,'primary');
        ii = ii +1;
    end

end

%% Prepare output
FRCinfo = struct();
FRCinfo.SSMorder   = order;
FRCinfo.SSMnonAuto = obj.Options.contribNonAuto;
FRCinfo.SSMispolar = ispolar;

% convert results to cell array
varargout{1} = FRC;

fdir = [pwd,'\data\',runid,'\SSMep.mat'];
save(fdir, 'FRC','FRCinfo');
end

function [] = continuation_BP(prim_runid,runid,lab,contOptions,cont_args,periodCell,ampCell, saveIC)
% Continues family of periodic orbits that emerges from a branch point
% setup coco
prob = coco_prob();
prob = cocoSet(contOptions, prob); %set default options 

prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', false);

% Continue from branch point
prob = ode_BP2po(prob, '', prim_runid, lab);

% store floquet multipliers in bd
prob = po_mult_add(prob,'po.orb'); 

%% Fix response period
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
[glueperiod,glueperiod_du] = deal(periodCell{:});

prob = coco_add_func(prob, 'OmegaT', glueperiod,glueperiod_du,  [],'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);


%% Monitor function for full system response
% Convert reduced results to full while continuing to avoid overhead of SSM
% computations

[ampData,ampNames,ampfunc] = deal(ampCell{:});
prob = coco_add_func(prob, 'amp', ampfunc, ampData, 'regular', ampNames,...
    'uidx', uidx([maps.x0_idx, maps.p_idx]));

% monitor function for saving the initial condition
if saveIC
    prob = coco_add_slot(prob, 'x0tobd', @init_cond, data, 'bddat');
end

%% Start continuation

fprintf('\n Run=''%s'': Continue secondary branch of periodic orbits in ''%s'' .\n', ...
    runid,prim_runid);

coco(prob, runid, [], cont_args{:});

end

function [data, y] = full_amplitude(~, data, u,obj)
% Saves the full steady state response 
outdof = data.outdof;
x0     = u(1:2);
state  = x0(1)+1j*x0(2);
eps    = u(end);

%% Reconstruct orbit
phi = linspace(0,2*pi,data.nt);
%reduced coordinate
orb = transpose(state).*exp(1i*phi);
red = [orb;conj(orb)];

%% Get full response
[z] = reduced_to_full(red,data.W,obj.W_1,eps);


% Save IC of full system in data
if data.saveIC
    data.IC = z(:,1);
end
%% Take the norm of the solution
% Frobenius
Znorm = norm(z,'fro')/sqrt(data.nt-1);

% at outdofs
Aout = [];
for k=1:numel(outdof)
    amp = norm(z(outdof(k),:),'inf');
    Aout = [Aout amp];
end

y = [Aout,Znorm].';
end

function [z0,t0] = get_initial_sol(z0,p0,odefun_in,nCycle,Omega,ispolar,periodsRatio)

T0 = periodsRatio*2*pi/Omega;
tf = nCycle*T0;

odefun = @(t,x) odefun_in (t,x,p0);
[~, z0_po] = ode15s(odefun , [0 tf], z0);                 % transient
options = odeset('RelTol', 1e-9, 'AbsTol',1e-9);
[t0, z0] = ode45(odefun, [0 T0], z0_po(end,:)', options); % steady state

%regularize solutions if polar
if ispolar
    error('polar not implemented')
end
end

function FRCy = array2structArray(FRCx)

om    = FRCx.om;
ep    = FRCx.ep;
st    = FRCx.st;
SNidx = FRCx.SNidx;
PDidx = FRCx.PDidx;
TRidx = FRCx.TRidx;
BPidx = FRCx.BPidx;
Aout  = FRCx.Aout_frc;
Zout  = FRCx.Zout_frc;
Znorm = FRCx.Znorm_frc;
Zic   = FRCx.Z0_frc;

saveIC = ~isempty(Zic);

numPts = numel(om);
FRCy   = cell(numPts,1);
for i=1:numPts
    FRC = struct( 'stability', st(i), 'Omega', om(i) ,...
    'epsilon', ep(i), 'Aout', Aout(i,:), 'Zout', [],...
    'Znorm', Znorm(i), 'Zic', [], 'isSN', false,'isPD',false,'isTR',false,'isBP',false);
    if saveIC
        FRC.Zic = Zic(i,:);
    end
    if ismember(i,SNidx)
        FRC.isSN = true;
    end
    
    if ismember(i,PDidx)
        FRC.isPD = true;
    end
    
    if ismember(i,TRidx)
        FRC.isTR = true;
    end
    
    if ismember(i,BPidx)
        FRC.isBP = true;
    end
    FRCy{i} = FRC;
end

FRCy = cat(1,FRCy{:});
end

function [data, res] = init_cond(prob, data, command, varargin) 
%init_cond   Append initial condition to BD.

res = {};
switch command
  case 'init'
    res = 'IC';
  case 'data'
    [fdata] = coco_get_func_data(prob, 'amp','data');
    res  = fdata.IC;
end
end
