coco_use_recipes_toolbox atlas1d_v5 % Add atlas1d_v5 atlas algorithm to search path

%% Example 14.3

% Test functionality that supports starting continuation from initial point
% on boundary of computational domain.

prob = coco_add_func(coco_prob(), 'circle', @circle, [], ...
  'zero', 'u0', [1.5; 1]);
prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);

coco(prob, '1', [], 1, {'x' 'y'}, [1 2]);     % Inside computational domain [Not in Recipes for Continuation]
coco(prob, '1', [], 1, {'x' 'y'}, [1.5 2]);   % Continue in admissible direction
coco(prob, '1', [], 1, {'x' 'y'}, [1.5 1.5]); % Terminate after locating initial chart

coco_use_recipes_toolbox % Remove atlas1d_v5 atlas algorithm from search path
