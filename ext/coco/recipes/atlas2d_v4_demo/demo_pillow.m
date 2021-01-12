coco_use_recipes_toolbox atlas2d_v4 % Add atlas2d_v4 atlas algorithm to search path

%% Example 13.2

% The continuation problem encoded below consists of a single zero function
% in terms of three continuation variables, a family of three monitor
% functions that evaluate to the continuation variables, and the
% corresponding inactive continuation parameters 'x', 'y', and 'z'. Its
% dimensional deficit equals -1. A two-dimensional surface is obtained by
% releasing 'x', 'y', and 'z' and allowing them to vary during
% continuation.

prob = coco_add_func(coco_prob(), 'pillow', @pillow, [], ...
  'zero', 'u0', [1; 0; 0]);
prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'NPR', 200);
prob = coco_set(prob, 'cont', 'h', 0.15, 'almax', 30);
coco(prob, 'pillow', [], 2, {'x' 'y' 'z'});

coco_use_recipes_toolbox % Add atlas2d_v4 atlas algorithm to search path
