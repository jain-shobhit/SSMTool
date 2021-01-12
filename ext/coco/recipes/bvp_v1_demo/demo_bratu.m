coco_use_recipes_toolbox coll_v1 bvp_v1 % add the coll_v1 and bvp_v1 toolboxes to the search path

%% Example 8.2 The Bratu boundary value problem

% The continuation problem encoded within the call to the coco entry-point
% function by the bvp_isol2seg constructor corresponds to a family of 101
% zero functions (80 collocation conditions, 18 continuity conditions, 3
% boundary conditions) in 102 continuation variables (100 basepoint values,
% 1 interval length, 1 problem parameter), a single monitor function, and
% the corresponding inactive continuation parameter 'p'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'p' is released and allowed to vary during continuation.

coll_args = {@brat, [0;1], zeros(2), 'p', 0};
bvp_args  = [coll_args, {@brat_bc, @brat_bc_DFDX}];
prob      = coco_prob();
coco(prob, 'brat1', @bvp_isol2seg, bvp_args{:}, 1, 'p', [0 4]);

% Use the following to graph a one-parameter bifurcation diagram for pairs
% of values of the problem parameter and the Euclidean norm of the
% vector of continuation variables.

% bd   = coco_bd_read('brat1');    % Extract bifurcation data cell array
% par  = coco_bd_col(bd, 'p');     % Extract column data
% nrmx = coco_bd_col(bd, '||U||'); % Extract column data
% plot(par, nrmx, 'b.-');

%% test restart of bvp [Not in Recipes for Continuation]

% The recursive parsing by the bvp_sol2seg constructor within the call to
% the coco entry-point function reconstructs the continuation problem
% considered previously using the solution and corresponding data
% associated with the 6th labeled solution from the previous run. The
% continuation parameter 'p' is again released and allowed to vary during
% continuation. The bifurcation data cell array is assigned to the bd
% variable.

prob = coco_prob();
bd = coco(prob, 'brat2', @bvp_sol2seg, 'brat1', 6, 1, 'p', [0 3.5]);

coco_use_recipes_toolbox % remove the coll_v1 and bvp_v1 toolboxes from the search path
