% minimize x s.t 2-x<=0

%% initial state in feasible region
prob = coco_prob();
prob = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'u0', 3);
prob = coco_add_func(prob, 'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', 1);
% adjoint
prob = coco_add_adjt(prob, 'bound', 'ncp.bound');
prob = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', 1);
coco(prob, 'run1', [], 1, {'obj' 'd.obj' 'ncp.bound' 'bound'}, [0 6]);

%% initial state in infeasible region
prob = coco_prob();
prob = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'u0', 1);
prob = coco_add_func(prob, 'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', 1);
% adjoint
prob = coco_add_adjt(prob, 'bound', 'ncp.bound');
prob = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', 1);
coco(prob, 'run2', [], 1, {'d.obj' 'obj' 'ncp.bound' 'bound'}, [0 1]);

%% restart from terminal point
bd   = coco_bd_read('run2');
labs = coco_bd_labs(bd, 'EP');

prob = coco_prob();
chart = coco_read_solution('bound', 'run2', labs(end), 'chart');
prob = coco_add_func(prob, 'bound', @bound, @bound_du, [], ...
  'inequality', 'bound', 'u0', chart.x);
prob = coco_add_func(prob, 'obj', @obj, @obj_du, [], ...
  'inactive', 'obj', 'uidx', 1);
% adjoint
chart = coco_read_adjoint('bound', 'run2', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'bound', 'ncp.bound', 'l0', chart.x);
chart = coco_read_adjoint('obj', 'run2', labs(end), 'chart');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', 1, 'l0', chart.x);
coco(prob, 'run3', [], 1, {'ncp.bound' 'obj' 'bound' 'd.obj'}, [0 2]);
