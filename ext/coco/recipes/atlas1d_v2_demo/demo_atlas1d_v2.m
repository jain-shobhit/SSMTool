coco_use_recipes_toolbox atlas1d_v2 % Add atlas1d_v2 atlas algorithm to search path

%% Example 12.2

% Check automatic step size refinement
prob = coco_add_func(coco_prob(), 'ellipse', @ellipse, [], ...
  'zero', 'u0', [0.7; 0.8]);
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'PtMX', 20, 'hmin', 0.065);
bd = coco(prob, 'run', [], 1);
coco_bd_col(bd, 'StepSize')

%% Further tests of conditional flows: theta [Not in Recipes for Continuation] 

prob = coco_prob();
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);

% theta == 0 : tangent predictor (explicit Euler)
pprob = coco_set(prob, 'cont', 'theta', 0, 'PtMX', 5);
pprob = coco_add_func(pprob, 'ellipse', @ellipse, [], 'zero', 'u0', [1.2;0.9] );
pprob = coco_add_pars(pprob, '', [1 2], {'x' 'y'});
coco(pprob, 'th1', [], 1, {'x' 'y'});

% theta = 0.5 : linearly implicit mid-point rule
pprob = coco_set(prob, 'cont', 'theta', 0.5, 'PtMX', 5);
pprob = coco_add_func(pprob, 'ellipse', @ellipse, [], 'zero', 'u0', [1.2;0.9] );
pprob = coco_add_pars(pprob, '', [1 2], {'x' 'y'});
coco(pprob, 'th2', [], 1, {'x' 'y'});

% theta = 1 : linearly implicit Euler
pprob = coco_set(prob, 'cont', 'theta', 1, 'PtMX', 5);
pprob = coco_add_func(pprob, 'ellipse', @ellipse, [], 'zero', 'u0', [1.2;0.9] );
pprob = coco_add_pars(pprob, '', [1 2], {'x' 'y'});
coco(pprob, 'th3', [], 1, {'x' 'y'});

% view log files of previous runs
coco_view_log('th1', 'type')
coco_view_log('th2', 'type')
coco_view_log('th3', 'type')

coco_use_recipes_toolbox % Remove atlas1d_v2 atlas algorithm from search path
