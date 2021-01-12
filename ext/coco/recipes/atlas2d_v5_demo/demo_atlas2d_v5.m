coco_use_recipes_toolbox atlas2d_v5 % Add atlas2d_v5 atlas algorithm to search path

%% Example 14.4

% The continuation problem structure encoded below corresponds to a single
% zero function in three continuation variables, a family of four monitor
% functions, and the corresponding active continuation parameter 'p' and
% inactive continuation parameters 'x', 'y', and 'z'. Its dimensional
% deficit equals -1. We obtain a two-dimensional solution manifold by releasing
% 'x', 'y', and 'z', and allowing the to vary during continuation. The
% computational domain boundary imposed by the limits on 'p' restrict
% continuation to a wedge in three-dimensional space.

prob = coco_add_func(coco_prob(), 'cylinder', @cylinder, [], ...
  'zero', 'u0', [1;0;0]+sqrt([0.5;0.55;0]) );
prob = coco_add_func(prob, 'angle', @angle, [], 'active', 'p', ...
  'uidx', (1:3)');
prob = coco_add_pars(prob, '', [1, 2, 3], {'x' 'y' 'z'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'h', 0.4, 'almax', 20, 'PtMX', 250);
coco(prob, 'cylinder', [], 2, {'x' 'z' 'y' 'p'}, ...
  {[], [], [], [-0.5 0.5]});

%% Sphere [Not in Recipes for Continuation]

% The continuation problem structure encoded below corresponds to a single
% zero function in three continuation variables, a family of three monitor
% functions that evaluate to the continuation variables, and the
% corresponding inactive continuation parameters 'x', 'y', and 'z'. Its
% dimensional deficit equals -1. We obtain a two-dimensional solution
% manifold by releasing 'x', 'y', and 'z' and allowing these to vary during
% continuation. The computational domain boundary imposed by the limits on
% 'x' restrict continuation to a layer in three-dimensional space.

prob = coco_add_func(coco_prob(), 'sphere', @sphere, [], 'zero', ...
  'u0', [1;1;0] );
prob = coco_add_pars(prob, '', 1:3, {'x', 'y', 'z'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'PtMX', 150, 'almax', 35, 'h', 0.8);
coco(prob, 'sphere', [], 2, {'x', 'y', 'z'}, {[0.5 1.5]});

% Use the following commands to visualize the two-dimensional manifold
% surface.
%
% atlas = coco_bd_read('sphere', 'atlas');
% plot_trisurf(atlas.charts, 1, 2, 3);

coco_use_recipes_toolbox % Remove atlas2d_v5 atlas algorithm to search path
