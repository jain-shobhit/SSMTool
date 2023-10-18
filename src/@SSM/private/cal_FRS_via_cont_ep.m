function cal_FRS_via_cont_ep(oid,fdata,obj,m,nCycle,scale_state,scale_obs,mFreqs,optData,outdof,optdof,parRange,varargin)
% This function calcualtes the forced response surface of the system via
% two-dimensional manifold continuation algorithms for fixed points of the
% reduced order models.

%% setup zero fucntions for fixed points

ispolar       = strcmp(obj.FRCOptions.coordinates, 'polar');
isbaseForce   = obj.System.Options.BaseExcitation;
initialSolver = obj.FRCOptions.initialSolver;
fdata.ispolar = ispolar;
fdata.isbaseForce = isbaseForce;
if ispolar
    odefun = @(z,p) ode_2mDSSM_polar(z,p,fdata);
    odefun_dfdx = @(z,p) ode_2mDSSM_polar_DFDX(z,p,fdata);
    odefun_dfdp = @(z,p) ode_2mDSSM_polar_DFDP(z,p,fdata);
else
    odefun = @(z,p) ode_2mDSSM_cartesian(z,p,fdata);
    odefun_dfdx = @(z,p) ode_2mDSSM_cartesian_DFDX(z,p,fdata);
    odefun_dfdp = @(z,p) ode_2mDSSM_cartesian_DFDP(z,p,fdata);    
end
funcs  = {odefun,odefun_dfdx,odefun_dfdp};

% 
%% setup continuation problem
coco_func_data.pointers('set', []); % only necessary for atlas_kd
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);
prob = cocoSet(obj.contOptions, prob);
% construct initial solution
if obj.System.order==2
    p0 = [obj.System.Omega; obj.System.fext.epsilon];
else
    p0 = [obj.System.Omega; obj.System.Fext.epsilon];
end
[p0,z0] = initial_fixed_point(p0,initialSolver,ispolar,...
    odefun,nCycle,m,varargin{:});

% call ep-toolbox
prob = ode_isol2ep(prob, '', funcs{:}, z0, {'om' 'eps'}, p0);
% define monitor functions to state variables
if isempty(scale_state)
    scales = ones(2*m,1);
else
    scales = scale_state;
end
[prob, args1, args2] = monitor_scaled_states(prob, ispolar, m, scales);
% observables
uidxpo  = optData.uidxpo;
nout    = numel(outdof);
argsoutLinf = cell(nout,1);
argsoutL2   = cell(nout,1);
outidx  = zeros(nout,1);
for k=1:nout
    argsoutLinf{k} = strcat('linfz',num2str(outdof(k)));
    argsoutL2{k}   = strcat('l2z',num2str(outdof(k)));
    idxk = find(outdof(k)==optdof);
    assert(numel(idxk)>0,'outdof is not included in optdof');
    outidx(k) = idxk;
end
optData.outdof = outdof;
optData.outidx = outidx;
optData.mFreqs = mFreqs;
if isempty(scale_obs)
    optData.scales = ones(2*nout+1,1);
else
    optData.scales = scale_obs;
end
prob = coco_add_func(prob, 'ampobj', @frs_output, optData, ...
    'active', [{'amp'},argsoutL2(:)',argsoutLinf(:)'], 'uidx', uidxpo);

runid = coco_get_id(oid, 'FRSep');

if ispolar
    cont_args = [{'om'}, {'eps'}, args1(:)', args2(:)', {'amp'}, argsoutL2(:)', argsoutLinf(:)'];
else
    % also monitor rhos
    args3 = cell(m,1);
    for k=1:m
        args3{k} = strcat('rho',num2str(k));
    end
    rhodata = struct();
    rhodata.scale = scales(1:2:end-1);
    prob = coco_add_func(prob, 'radius', @cal_rhos,...
        rhodata, 'active', args3(:)', 'uidx', 1:2*m);
    cont_args = [{'om'}, {'eps'}, args1(:)', args2(:)', {'amp'}, argsoutL2(:)', argsoutLinf(:)', args3(:)'];
end

fprintf('\n Run=''%s'': Continue equilibria along FRS.\n', ...
  runid);

coco(prob, runid, [], 2, cont_args, parRange);

figure; grid on
atlas = coco_bd_read(runid, 'atlas');
plot_atlas_kd(atlas.charts,1,2,3);
view([20,50]);

end