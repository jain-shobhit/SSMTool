coco_use_recipes_toolbox atlas2d_v2 % Add atlas2d_v2 atlas algorithm to search path

%% Example 13.1

% The continuation problem structure encoded below corresponds to a single
% zero function in terms of three continuation variables, a family of three
% monitor functions that evaluate to the continuation variables, and the
% corresponding inactive continuation parameters 'x', 'y', and 'z'. Its
% dimensional deficit is -1. A two-dimensional solution manifold results by
% releasing 'x', 'y', and 'z' and allowing these to vary during
% continuation.

% In the first run, the number of available directions for
% continuation from individual charts is 6, but the algorithm terminates
% prematurely after going once around a great circle. In the second run,
% the number of available directions for continuation is only 4. The
% algorithm covers the entire sphere without redundancy.

prob = coco_add_func(coco_prob(), 'sphere', @sphere, [], ...
  'zero', 'u0', [2;0;0]);
prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'h', 0.5, 'almax', 35); % Number of available directions = 6
coco(prob, 'sphere1', [], 2, {'x' 'y' 'z'});

prob = coco_set(prob, 'cont', 'Ndirs', 4, 'PtMX', 200); % Number of available directions = 4
coco(prob, 'sphere2', [], 2, {'x' 'y' 'z'});

coco_use_recipes_toolbox % Remove atlas2d_v2 atlas algorithm from search path
