coco_use_recipes_toolbox alg_v5 % add the alg_v5 toolbox to the search path

%% Example 4.7

% The continuation problem encoded below corresponds to a family of two
% zero function in terms of four continuation variables, a family of two
% monitor functions that evaluates to the second and fourth continuation
% variables, respectively, and two corresponding inactive continuation
% parameters 'y1' and 'y2'. Its dimensional deficit equals 0. The call to
% the coco entry-point function indicates a desired manifold dimension of
% 1. To this end, the continuation parameter 'y1' is released and allowed
% to vary during continuation. The continuation parameter 'y2' remains
% inactive during continuation.

prob = coco_prob();
fhan = @(x,p) x^2+(p-1)^2-1;
prob = alg_isol2eqn(prob, 'eqn1', fhan, 0.9, 'y1', 1.1);
prob = alg_isol2eqn(prob, 'eqn2', fhan, 0.1, 'y2', 1.9);
coco(prob, 'run1', [], 1, {'y1' 'y2'}, [0.5 2]);

alg_read_solution('eqn1', 'run1', 4)
alg_read_solution('eqn2', 'run1', 4)

%% Example 4.8

% The continuation problem encoded below is identical to that above, but
% constructed from a stored solution from the previous run.

prob = coco_prob();
prob = alg_sol2eqn(prob, 'eqn3', 'run1', 'eqn1', 4);
prob = alg_sol2eqn(prob, 'eqn4', 'run1', 'eqn2', 4);
coco(prob, 'run2', [], 1, {'y2' 'y1'}, [0.5 2]);

coco_use_recipes_toolbox % remove the alg_v5 toolbox from the search path
