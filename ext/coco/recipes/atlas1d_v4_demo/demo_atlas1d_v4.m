coco_use_recipes_toolbox atlas1d_v4 % Add atlas1d_v4 atlas algorithm to search path

%% Example 14.2

% Continuation along both directions of the solution manifold
prob = coco_add_func(coco_prob(), 'circle', @circle, [], ...
  'zero', 'u0', [1.5; 1]);
prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
coco(prob, '1', [], 1, {'x' 'y'}, [1 2]);

%% Further tests of functionality: initial point on domain boundary [Not in Recipes for Continuation]
coco(prob, '1', [], 1, {'x' 'y'}, [1.5 2]);   % Available direction into computational domain
coco(prob, '1', [], 1, {'x' 'y'}, [1.5 1.5]); % No available direction into computational domain

coco_use_recipes_toolbox % Remove atlas1d_v4 atlas algorithm from search path
