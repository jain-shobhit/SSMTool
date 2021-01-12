coco_use_recipes_toolbox alg_v4 % add the alg_v4 toolbox to the search path

%% Example 4.4

% The continuation problem encoded below corresponds to a single zero
% function in terms of two continuation variables, a single monitor
% function that evaluates to the second continuation variables, and a
% corresponding inactive continuation parameter 'y'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'y' is released and allowed to vary during continuation.

prob = alg_isol2eqn(@(x,p) x^2+(p-1)^2-1, 0.9, 'y', 1.1);
coco(prob, 'run1', [], 1, 'y', [0.1 5]);

%% Example 4.5

% The continuation problem encoded below is identical to that above, but
% constructed from a stored solution from the previous run.

prob = alg_sol2eqn('run1', 9);
coco(prob, 'run2', [], 1, 'y', [0.1 1]);

sol = alg_read_solution('run2', 2)

coco_use_recipes_toolbox % remove the alg_v4 toolbox from the search path
