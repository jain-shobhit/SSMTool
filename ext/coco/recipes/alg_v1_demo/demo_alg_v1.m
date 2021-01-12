coco_use_recipes_toolbox alg_v1 % add the alg_v1 toolbox to the search path

%% Example 4.1

% The continuation problem encoded below corresponds to a single zero
% function in terms of two continuation variables, a single monitor
% function that evaluates to the second continuation variables, and a
% corresponding inactive continuation parameter 'y'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'y' is released and allowed to vary during continuation.

prob = alg_construct_eqn(@circle, 1, 'y', 1.1);
coco(prob, 'run1', [], 1, 'y', [0.1 5]);

%% Example 4.1 with anonymous function

% The continuation problem encoded below is identical to that above, other
% than the use of the anonymous function as the first input argument to the
% alg_construct_eqn toolbox constructor.

prob = alg_construct_eqn(@(x,p) x^2+(p-1)^2-1, 1, 'y', 1.1);
coco(prob, 'run1', [], 1, 'y', [0.1 5]);

coco_use_recipes_toolbox % remove the alg_v1 toolbox from the search path
