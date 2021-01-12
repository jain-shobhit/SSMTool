% minimize x^2+y s.t 2-y<=0

%% initial state in feasible region
prob = coco_prob();
prob = coco_add_pars(prob, 'x', 'x', 1);
prob = coco_add_pars(prob, 'y', 'y', 3);
prob = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'uidx', 2);
prob = coco_add_func(prob, 'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', [1 2]);
% adjoint
prob = coco_add_adjt(prob, 'x', 'd.x');
prob = coco_add_adjt(prob, 'y', 'd.y');
prob = coco_add_adjt(prob, 'bound', 'ncp.bound', 'aidx', 2);
prob = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1 2]);

cont_pars = {'obj' 'x' 'd.obj' 'd.y' 'y' 'bound' 'd.x' 'ncp.bound'};
coco(prob, 'run1', [], 1, cont_pars, [0 6]);

%% branch switch
bd  = coco_bd_read('run1');
lab = coco_bd_labs(bd, 'BP');

prob  = coco_prob();
chart = coco_read_solution('x', 'run1', lab, 'chart');
prob  = coco_add_pars(prob,'x', 'x', chart.x);
chart = coco_read_solution('y', 'run1', lab, 'chart');
prob  = coco_add_pars(prob,'y', 'y', chart.x);
prob  = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'uidx', 2);
prob  = coco_add_func(prob, 'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', [1 2]);
% adjoint
cdata = coco_get_chart_data(chart, 'lsol');

[chart, lidx] = coco_read_adjoint('x', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'x', 'd.x', 'l0', chart.x, ...
  'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('y', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'y', 'd.y', 'l0', chart.x, ...
  'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('bound', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'bound', 'ncp.bound', 'aidx', 2, ...
  'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('obj', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1 2], ...
  'l0', chart.x, 'tl0', cdata.v(lidx));

cont_pars = {'d.obj' 'obj' 'x' 'd.y' 'y' 'bound' 'd.x' 'ncp.bound'};
coco(prob, 'run2', [], 1, cont_pars, [0 1]);

%% drive d.y to zero
bd   = coco_bd_read('run2');
labs = coco_bd_labs(bd, 'EP');

prob  = coco_prob();
chart = coco_read_solution('x', 'run2', labs(end), 'chart');
prob  = coco_add_pars(prob,'x', 'x', chart.x);
chart = coco_read_solution('y', 'run2', labs(end), 'chart');
prob  = coco_add_pars(prob,'y', 'y', chart.x);
prob  = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'uidx', 2);
prob  = coco_add_func(prob, 'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', [1 2]);
% adjoint
chart = coco_read_adjoint('x', 'run2', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'x', 'd.x', 'l0', chart.x);
chart = coco_read_adjoint('y', 'run2', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'y', 'd.y', 'l0', chart.x);
chart = coco_read_adjoint('bound', 'run2', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'bound', 'ncp.bound', 'aidx', 2, ...
  'l0', chart.x);
chart = coco_read_adjoint('obj', 'run2', labs(end), 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1 2], 'l0', chart.x);

cont_pars = {'d.y' 'obj' 'x'  'ncp.bound' 'y' 'bound' 'd.x' 'd.obj'};
coco(prob, 'run3', [], 1, cont_pars, [-1 0]);

%% drive ncp.bound to 0
bd   = coco_bd_read('run3');
labs = coco_bd_labs(bd, 'EP');

prob  = coco_prob();
chart = coco_read_solution('x', 'run3', labs(end), 'chart');
prob  = coco_add_pars(prob,'x', 'x', chart.x);
chart = coco_read_solution('y', 'run3', labs(end), 'chart');
prob  = coco_add_pars(prob,'y', 'y', chart.x);
prob  = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'uidx', 2);
prob  = coco_add_func(prob,'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', [1 2]);
% adjoint
chart = coco_read_adjoint('x', 'run3', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'x', 'd.x', 'l0', chart.x);
chart = coco_read_adjoint('y', 'run3', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'y', 'd.y', 'l0', chart.x);
chart = coco_read_adjoint('bound', 'run3', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'bound', 'ncp.bound', 'aidx', 2, ...
  'l0', chart.x);
chart = coco_read_adjoint('obj', 'run3', labs(end), 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1 2], 'l0', chart.x);

cont_pars = {'ncp.bound' 'obj' 'x'  'y' 'bound' 'd.x' 'd.y' 'd.obj'};
coco(prob, 'run3', [], 1, cont_pars, [sqrt(2)-2 0]);

%% initial state in infeasible region
prob = coco_prob();
prob = coco_add_pars(prob, 'x', 'x', 1);
prob = coco_add_pars(prob, 'y', 'y', 1);
prob = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'uidx', 2);
prob = coco_add_func(prob,'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', [1 2]);
% adjoint
prob = coco_add_adjt(prob, 'x', 'd.x');
prob = coco_add_adjt(prob, 'y', 'd.y');
prob = coco_add_adjt(prob, 'bound', 'ncp.bound', 'aidx', 2);
prob = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1 2]);

cont_pars = {'obj' 'x' 'd.obj' 'd.y' 'y' 'bound' 'd.x' 'ncp.bound'};
coco(prob, 'run1', [], 1, cont_pars, [0 6]);

%% branch switch
bd  = coco_bd_read('run1');
lab = coco_bd_labs(bd, 'BP');

prob  = coco_prob();
chart = coco_read_solution('x', 'run1', lab, 'chart');
prob  = coco_add_pars(prob,'x', 'x', chart.x);
chart = coco_read_solution('y', 'run1', lab, 'chart');
prob  = coco_add_pars(prob,'y', 'y', chart.x);
prob  = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'uidx', 2);
prob  = coco_add_func(prob, 'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', [1 2]);
% adjoint
cdata = coco_get_chart_data(chart, 'lsol');

[chart, lidx] = coco_read_adjoint('x', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'x', 'd.x', 'l0', chart.x, ...
  'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('y', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'y', 'd.y', 'l0', chart.x, ...
  'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('bound', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'bound', 'ncp.bound', 'aidx', 2, ...
  'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('obj', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1 2], ...
  'l0', chart.x, 'tl0', cdata.v(lidx));

cont_pars = {'d.obj' 'obj' 'x' 'd.y' 'y' 'bound' 'd.x' 'ncp.bound'};
coco(prob, 'run2', [], 1, cont_pars, [0 1]);

%% drive d.y to 0
bd   = coco_bd_read('run2');
labs = coco_bd_labs(bd, 'EP');

prob  = coco_prob();
chart = coco_read_solution('x', 'run2', labs(end), 'chart');
prob  = coco_add_pars(prob,'x', 'x', chart.x);
chart = coco_read_solution('y', 'run2', labs(end), 'chart');
prob  = coco_add_pars(prob,'y', 'y', chart.x);
prob  = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'uidx', 2);
prob  = coco_add_func(prob, 'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', [1 2]);
% adjoint
chart = coco_read_adjoint('x', 'run2', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'x', 'd.x', 'l0', chart.x);
chart = coco_read_adjoint('y', 'run2', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'y', 'd.y', 'l0', chart.x);
chart = coco_read_adjoint('bound', 'run2', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'bound', 'ncp.bound', 'aidx', 2, ...
  'l0', chart.x);
chart = coco_read_adjoint('obj', 'run2', labs(end), 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1 2], 'l0', chart.x);

cont_pars = {'d.y' 'obj' 'x'  'y' 'bound' 'd.x' 'ncp.bound' 'd.obj'};
coco(prob, 'run3', [], 1, cont_pars, [-1 0]);

%% drive ncp.bound to zero
bd  = coco_bd_read('run3');
lab = coco_bd_labs(bd, 'EP');

prob  = coco_prob();
chart = coco_read_solution('x', 'run3', labs(end), 'chart');
prob  = coco_add_pars(prob,'x', 'x', chart.x);
chart = coco_read_solution('y', 'run3', labs(end), 'chart');
prob  = coco_add_pars(prob,'y', 'y', chart.x);
prob  = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'uidx', 2);
prob  = coco_add_func(prob, 'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', [1 2]);
% adjoint
chart = coco_read_adjoint('x', 'run3', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'x', 'd.x', 'l0', chart.x);
chart = coco_read_adjoint('y', 'run3', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'y', 'd.y', 'l0', chart.x);
chart = coco_read_adjoint('bound', 'run3', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'bound', 'ncp.bound', 'aidx', 2, ...
'l0', chart.x);
chart = coco_read_adjoint('obj', 'run3', labs(end), 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1 2], 'l0', chart.x);

cont_pars = {'ncp.bound' 'obj' 'x'  'y' 'bound' 'd.x' 'd.y' 'd.obj'};
coco(prob, 'run4', [], 1, cont_pars, [0 2]);
