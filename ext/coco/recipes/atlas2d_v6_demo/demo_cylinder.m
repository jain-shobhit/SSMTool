coco_use_recipes_toolbox atlas2d_v6 % Add atlas2d_v6 atlas algorithm to search path

%% Example 14.5

% The continuation problem structure encoded below corresponds to a single
% zero function in three continuation variables, a family of three monitor
% functions that evaluate to the continuation variables, and the
% corresponding inactive continuation parameters 'x', 'y', and 'z'. Its
% dimensional deficit equals -1. We obtain a two-dimensional solution
% manifold by releasing 'x', 'y', and 'z', and allowing the to vary during
% continuation. In the first run, the initial point lies on a boundary of
% the computation domain imposed by a single interval limit on 'z'. In the
% second run, the initial point lies on the intersection of two boundaries
% of the computational domain imposed by an interval limit on 'x' and an
% interval limit on 'z', respectively.

prob = coco_add_func(coco_prob(), 'cylinder', @cylinder, [], ...
  'zero', 'u0', [1; -1; 0]);
prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'Ndirs', 10, 'PtMX', 25);
coco(prob, 'cylinder1', [], 2, {'x' 'y' 'z'}, {[] [] [0 1]});
coco(prob, 'cylinder2', [], 2, {'x' 'y' 'z'}, {[0 1] [] [0 1]});

%% Additional tests [Not in Recipes for Continuation]

% The continuation problem encoded below is identical to that above. In
% this run, the computational domain is bounded by interval limits on all
% three continuation parameters, and the initial point lies on an
% intersection of two such boundaries.

prob = coco_add_func(coco_prob(), 'cylinder', @cylinder, [], ...
  'zero', 'u0', [1; -1; -1]);
prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'Ndirs', 10, 'PtMX', 200, 'h', .2);
coco(prob, 'cylinder', [], 2, {'y' 'z' 'x'}, ...
  {[-2 0] [-1 0] [-1 1]});

% The continuation problem encoded below is identical to that above. In
% this run, the computational domain is bounded by interval limits on 'z'
% and the initial point lies outside the computational domain

prob = coco_add_func(coco_prob(), 'cylinder', @cylinder, [], 'zero', ...
  'x0', [2;0;-2] );
prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'h', .5, 'almax', 35, 'PtMX', 200);
coco(prob, 'cylinder', [], 2, {'y' 'z' 'x'}, {[] [-1 1] []});

coco_use_recipes_toolbox % Remove atlas2d_v6 atlas algorithm from search path
