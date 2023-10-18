function [runidom,runidbc,args1,args2] = damped_backbone_linf(oid,contOptions,Options,funcs,z0,p0,ispolar,optData,parRange)
% DAMPED_BACKBONE_LINF This function extracts ridges and trenches of a FRS
% whose amplitude is defined in terms of LINF norm. These ridges and trenches
% are defined as damped backbone curve in an unified way


prob = coco_prob();
prob = cocoSet(contOptions, prob);
% call ep-toolbox
prob = ode_isol2ep(prob, '', funcs{:}, z0, {'om' 'eps'}, p0);
% define monitor functions to state variables
m = numel(z0)/2;
[prob, args1, args2] = monitor_states(prob, ispolar, m);
% define optimization objective
objfun = {@po_amp @po_amp_du};
uidxpo = optData.uidxpo;
prob = coco_add_func(prob, 'optobj', objfun{:}, optData, ...
    'inactive', 'obj', 'uidx', uidxpo, 'u0', 0);
pidx = coco_get_func_data(prob, 'optobj', 'uidx');
pidx = pidx(end);
prob = coco_add_pars(prob, 'time', pidx, 't');
% adjoints
prob = adjt_isol2ep(prob, ''); % adjoint
prob = coco_add_adjt(prob, 'optobj', 'd.obj', 'aidx', uidxpo);
aidx = coco_get_adjt_data(prob, 'optobj', 'axidx');
prob = coco_add_adjt(prob, 'time', 'd.t', 'aidx', aidx(end));
runidt = coco_get_id(oid, 'contPhase');

fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  runidt);
cont_args = [{'obj' 'd.obj' 't' 'd.om' 'd.eps'},args1(:)',args2(:)',{'om'},{'eps'}]; 
bd0 = coco(prob, runidt, [], 1, cont_args, {[],[],[0 2*pi/p0(1)]});

%% branch switching to follow the secondary branch
BPlab = coco_bd_labs(bd0, 'BP');
% find the one with maximum manginitude
objv  = coco_bd_col(bd0, 'obj');
BPidx = coco_bd_idxs(bd0, 'BP');
[~,idx] = max(abs(objv(BPidx)));
BPlab = BPlab(idx);
steps = Options.DBCstepFactor;

% zero problem
prob = coco_prob();
prob = cocoSet(contOptions, prob);
prob = coco_set(prob, 'cont', 'h_max', steps(1)*contOptions.h_max);
prob = ode_BP2ep(prob, '', runidt, BPlab);
[prob, ~, ~] = monitor_states(prob, ispolar, m);
chart = coco_read_solution('optobj',runidt, BPlab, 'chart');
prob = coco_add_func(prob, 'optobj', objfun{:}, optData, ...
    'inactive', 'obj', 'uidx', uidxpo, 'u0', chart.x(end));
pidx = coco_get_func_data(prob, 'optobj', 'uidx');
pidx = pidx(end);
prob = coco_add_pars(prob, 'time', pidx, 't');

% branch switch data
chart = coco_read_solution(runidt, BPlab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');
% adjoint
prob = adjt_BP2ep(prob, '', runidt, BPlab);
[chart, lidx] = coco_read_adjoint('optobj', runidt, BPlab, ...
  'chart', 'lidx');  
prob = coco_add_adjt(prob, 'optobj', 'd.obj', 'aidx', uidxpo, ...
    'l0', chart.x, 'tl0', cdata.v(lidx));
aidx = coco_get_adjt_data(prob, 'optobj', 'axidx');
[chart, lidx] = coco_read_adjoint('time', runidt, BPlab, ...
  'chart', 'lidx');  
prob = coco_add_adjt(prob, 'time', 'd.t', 'aidx', aidx(end), ...
    'l0', chart.x, 'tl0', cdata.v(lidx));

runidBP = coco_get_id(oid, 'timeBP');

% continuation
fprintf('\n Run=''%s'': Continue until d.obj=1.\n', ...
    runidBP);
temp = cont_args{1}; cont_args{1} = cont_args{2}; cont_args{2} = temp;
bd1 = coco(prob, runidBP, [], 1, cont_args,[0,1]);

%% continuation in om with d.obj=1 
EPlab = coco_bd_labs(bd1,'EP');
EPlab = max(EPlab);
% zero problem
prob = coco_prob();
prob = cocoSet(contOptions, prob);
prob = coco_set(prob, 'cont', 'h_max', steps(2)*contOptions.h_max);
prob = ode_ep2ep(prob, '', runidBP, EPlab);
[prob, ~, ~] = monitor_states(prob, ispolar, m);
chart = coco_read_solution('optobj',runidBP, EPlab, 'chart');
prob = coco_add_func(prob, 'optobj', objfun{:}, optData, ...
    'inactive', 'obj', 'uidx', uidxpo, 'u0', chart.x(end));
pidx = coco_get_func_data(prob, 'optobj', 'uidx');
pidx = pidx(end);
prob = coco_add_pars(prob, 'time', pidx, 't');

% adjoint
prob = adjt_ep2ep(prob, '', runidBP, EPlab);
chart = coco_read_adjoint('optobj', runidBP, EPlab, 'chart');  
prob = coco_add_adjt(prob, 'optobj', 'd.obj', 'aidx', uidxpo, 'l0', chart.x);
aidx = coco_get_adjt_data(prob, 'optobj', 'axidx');
chart = coco_read_adjoint('time', runidBP, EPlab, 'chart');  
prob = coco_add_adjt(prob, 'time', 'd.t', 'aidx', aidx(end), 'l0', chart.x);

prob = coco_add_event(prob, 'OPT', 'd.om', '=', 0);

runidom = coco_get_id(oid, 'omega');
% continuation
fprintf('\n Run=''%s'': Continue in (t,omega) with d.obj=1.\n', ...
    runidom);
temp = cont_args{1};
cont_args{1} = cont_args{4}; cont_args{4} = cont_args{end-1}; cont_args{end-1} = temp;
bd2 = coco(prob, runidom, [], 1, cont_args,{[],[],[],parRange{1}});

%% continuation in (epf,om) with d.obj=1 and d.om=0    
OPTlab = coco_bd_labs(bd2, 'OPT');
numOPT = numel(OPTlab);
runidbc = cell(numOPT,1);
temp = cont_args{1}; cont_args{1} = cont_args{5};
cont_args{5} = cont_args{end}; cont_args{end} = temp;
for idbp = 1:numOPT
    disp(['Backbone curve along optimal point (OPT) ', num2str(OPTlab(idbp))])
    % zero problem
    prob = coco_prob();
    prob = cocoSet(contOptions, prob);
    prob = coco_set(prob, 'cont', 'h_max', steps(2)*contOptions.h_max);
    prob = ode_ep2ep(prob, '', runidom, OPTlab(idbp));
    [prob, ~, ~] = monitor_states(prob, ispolar, m);
    chart = coco_read_solution('optobj',runidom, OPTlab(idbp), 'chart');
    prob = coco_add_func(prob, 'optobj', objfun{:}, optData, ...
        'inactive', 'obj', 'uidx', uidxpo, 'u0', chart.x(end));
    pidx = coco_get_func_data(prob, 'optobj', 'uidx');
    pidx = pidx(end);
    prob = coco_add_pars(prob, 'time', pidx, 't');

    % adjoint
    prob = adjt_ep2ep(prob, '', runidom, OPTlab(idbp));
    chart = coco_read_adjoint('optobj', runidom, OPTlab(idbp), 'chart');  
    prob = coco_add_adjt(prob, 'optobj', 'd.obj', 'aidx', uidxpo, 'l0', chart.x);
    aidx = coco_get_adjt_data(prob, 'optobj', 'axidx');
    chart = coco_read_adjoint('time', runidom, OPTlab(idbp), 'chart');  
    prob = coco_add_adjt(prob, 'time', 'd.t', 'aidx', aidx(end), 'l0', chart.x);

    runidbc{idbp} = coco_get_id(oid, ['freq.bc',num2str(idbp)]);
    
%     prob = coco_add_event(prob, 'BC', 'BP','eps', '<', 0);

    % continuation
    fprintf('\n Run=''%s'': Continue in (t,omega, epf) with d.obj=1 and d.om=0.\n', ...
        runidbc{idbp});
    coco(prob, runidbc{idbp}, [], 1, cont_args,{[],[],[],parRange{1},parRange{2}});
end


end