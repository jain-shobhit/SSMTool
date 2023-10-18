function [runid,runidbc,args1,args2] = damped_backbone_l2(oid,contOptions,Options,funcs,z0,p0,ispolar,optData,parRange)
% DAMPED_BACKBONE_L2 This function extracts ridges and trenches of a FRS
% whose amplitude is defined in terms of L2 norm. These ridges and trenches
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
    'inactive', 'obj', 'uidx', uidxpo);
% adjoints
prob = adjt_isol2ep(prob, ''); % adjoint
prob = coco_add_adjt(prob, 'optobj', 'd.obj', 'aidx', uidxpo);

runid = coco_get_id(oid, 'freq.frc');
fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  runid);
cont_args = [{'obj' 'd.obj' 'om' 'd.eps'}, args1(:)' ,args2(:)',{'eps'}];
bd1 = coco(prob, runid, [], 1, cont_args, {[],[],parRange{1}});

%% branch switching to follow the secondary branch
BPlab = coco_bd_labs(bd1, 'BP');
steps = Options.DBCstepFactor;

%% for-loop over branch points
numBP = numel(BPlab);
runidbp = cell(numBP,1);
runidbc = cell(numBP,1);
for idbp = 1:numBP
    disp(['Backbone curve along branch point (BP) ', num2str(BPlab(idbp))])
    % branch switch data
    chart = coco_read_solution(runid, BPlab(idbp), 'chart');
    cdata = coco_get_chart_data(chart, 'lsol');

    % zero problem
    prob = coco_prob();
    prob = cocoSet(contOptions, prob);
    prob = coco_set(prob, 'cont', 'h_max', steps(1)*contOptions.h_max);
    prob = ode_BP2ep(prob, '', runid, BPlab(idbp));
    [prob, ~, ~] = monitor_states(prob, ispolar, m);
    cont_args = [{'d.obj' 'obj' 'om' 'd.eps'}, args1(:)' ,args2(:)',{'eps'}];
    prob = coco_add_func(prob, 'optobj', objfun{:}, optData, 'inactive', 'obj', 'uidx', uidxpo);

    % adjoint
    prob = adjt_BP2ep(prob, '', runid, BPlab(idbp));
    [chart, lidx] = coco_read_adjoint('optobj', runid, BPlab(idbp), ...
      'chart', 'lidx');  
    prob = coco_add_adjt(prob, 'optobj', 'd.obj', 'aidx', uidxpo, ...
      'l0', chart.x, 'tl0', cdata.v(lidx));

    runidbp{idbp} = coco_get_id(oid, ['freq.bp',num2str(idbp)]);

    % continuation
    fprintf('\n Run=''%s'': Continue until d.obj=1.\n', ...
        runidbp{idbp});
    bd2 = coco(prob, runidbp{idbp}, [], 1, cont_args,[0,1]);

    %% continuation along backbone curve by releasing epsilon
    EPlab = coco_bd_labs(bd2, 'EP');
    EPlab = max(EPlab);
    % zero problem
    prob = coco_prob();
    prob = cocoSet(contOptions, prob);
    prob = coco_set(prob, 'cont', 'h_max', steps(2)*contOptions.h_max);
    if ~isempty(Options.PtMXBCrun)
        prob = coco_set(prob, 'cont', 'PtMX', Options.PtMXBCrun);
    end
    prob = ode_ep2ep(prob, '', runidbp{idbp}, EPlab);
    [prob, ~, ~] = monitor_states(prob, ispolar, m);
    cont_args = [{'d.eps' 'obj' 'om' 'eps'}, args1(:)' ,args2(:)',{'d.obj'}];
    prob = coco_add_func(prob, 'optobj', objfun{:}, optData, 'inactive', 'obj', 'uidx', uidxpo);

    % adjoint
    prob  = adjt_ep2ep(prob, '', runidbp{idbp}, EPlab);
    chart = coco_read_adjoint('optobj', runidbp{idbp}, EPlab, 'chart');  
    prob  = coco_add_adjt(prob, 'optobj', 'd.obj', 'aidx', uidxpo, ...
      'l0', chart.x);

    runidbc{idbp} = coco_get_id(oid, ['freq.bc',num2str(idbp)]);

    % continuation
    fprintf('\n Run=''%s'': Continue along backbone curve.\n', ...
        runidbc{idbp});    
    coco(prob, runidbc{idbp}, [], 1, cont_args,{[],[],parRange{1},parRange{2}});
end


end