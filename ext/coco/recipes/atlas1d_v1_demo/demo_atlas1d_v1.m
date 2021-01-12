coco_use_recipes_toolbox atlas1d_v1 % Add atlas1d_v1 atlas algorithm to search path

%% Section 12.1.2

% Each run below provides a test of functionality and flow of the
% atlas1d_v1 atlas algorithm.

% Regular points
prob = coco_add_func(coco_prob(), 'circle', @circle, [], ...
  'zero', 'u0', [1.5; 1]);
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'PtMX', 30);
coco(prob, 'run', [], 1);

prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
coco(prob, 'run', [], 1, {'x' 'y'});

% Failure to converge during continuation. [Typo in Eqn. (12.11) in Recipes
% for Continuation, 1st edition]
data.MX = true;
prob = coco_add_func(coco_prob(), 'func', @circle, data, ...
  'zero', 'u0', [1.5; 1]);
prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
coco(prob, 'run', [], 1, {'x' 'y'});

% Failure to converge on initial point.
prob = coco_add_func(coco_prob(), 'func', @empty, [], ...
  'zero', 'u0', [1; 1]);
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
coco(prob, 'run', [], 1);

% Singular Jacobian
prob = coco_add_func(coco_prob(), 'func', @singular, [], ...
  'zero', 'u0', [1; 1]);
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
coco(prob, 'run', [], 1);

% Starting from vertical tangent [Not in Recipes for Continuation]
prob = coco_add_func(coco_prob(), 'func', @circle, [], ...
  'zero', 'u0', [2; 0]);
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
coco(prob, 'run', [], 1);

% Starting from singular point
prob = coco_add_func(coco_prob(), 'lemniscate', @lemniscate, ...
  [], 'zero', 'u0', [0; 0]);
prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
coco(prob, 'run', [], 1, {'x' 'y'});

% Starting from singular point with tangent vector. No output from NWTN as
% first point is accepted without correction.
prob = coco_add_func(coco_prob(), 'lemniscate', @lemniscate, ...
  [], 'zero', 'u0', [0; 0], 't0', [1; 1]);
prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
coco(prob, 'run', [], 1, {'x' 'y'});


%% Further tests of conditional flows: PtMX [Not in Recipes for Continuation]

prob = coco_prob();
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);

% PtMX >= 0
pprob = coco_set(prob, 'cont', 'PtMX', 30);
pprob = coco_add_func(pprob, 'circle', @circle, [], 'zero', 'u0', [1.5;1] );
pprob = coco_add_pars(pprob, '', [1 2], {'x' 'y'});
coco(pprob, 'run', [], 1, {'x' 'y'});
% PtMX == 0
pprob = coco_set(prob, 'cont', 'PtMX', 0);
pprob = coco_add_func(pprob, 'circle', @circle, [], 'zero', 'u0', [1.5;1] );
pprob = coco_add_pars(pprob, '', [1 2], {'x' 'y'});
coco(pprob, 'run', [], 1, {'x' 'y'});
% PtMX <= 0
pprob = coco_set(prob, 'cont', 'PtMX', -30);
pprob = coco_add_func(pprob, 'circle', @circle, [], 'zero', 'u0', [1.5;1] );
pprob = coco_add_pars(pprob, '', [1 2], {'x' 'y'});
coco(pprob, 'run', [], 1, {'x' 'y'});

coco_use_recipes_toolbox % Remove atlas1d_v1 atlas algorithm from search path
