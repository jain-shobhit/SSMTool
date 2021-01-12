coco_use_recipes_toolbox atlas1d_v3 % Add atlas1d_v3 atlas algorithm to search path

%% Section 12.3.2

% Each run below provides a test of functionality and flow of the
% atlas1d_v3 atlas algorithm.

% Terminates when boundary becomes empty
prob = coco_add_func(coco_prob(), 'circle', @circle, [], ...
  'zero', 'u0', [1.5; 1]);
prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'PtMX', 100);
coco(prob, 'run', [], 1, {'x' 'y'});

% Terminates when gap between successive charts
prob = coco_add_func(coco_prob(), 'ellipse', @ellipse, [], ...
  'zero', 'u0', [0.6; 0.5]);
prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
coco(prob, 'run', [], 1, {'x' 'y'});

%% Example 14.1

% Specified computational domain
prob = coco_add_func(coco_prob(), 'circle', @circle, [], ...
  'zero', 'u0', [1.5; 1]);
prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
coco(prob, 'run', [], 1, {'x' 'y'}, [2 3]); % Outside computational domain
coco(prob, 'run', [], 1, {'x' 'y'}, [1 2]); % Inside computationl domain

coco_use_recipes_toolbox % Remove atlas1d_v1 atlas algorithm from search path
