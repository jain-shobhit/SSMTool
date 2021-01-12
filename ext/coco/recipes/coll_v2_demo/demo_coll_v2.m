coco_use_recipes_toolbox coll_v2 bvp_v1 % add the coll_v2 and bvp_v1 toolboxes to the search path

%% Continuation from zero-length segment [Not in Recipes for Continuation]

% The continuation problem structure encoded below corresponds to a family
% of 150 zero functions (120 collocation conditions, 27 continuity
% conditions, and three boundary conditions) in terms of 152 continuation
% variables (150 basepoint values, one interval length, and one problem
% parameter), a family of two monitor functions that evaluate to the
% problem parameter and the interval length, respectively, and the
% corresponding inactive continuation parameters 'p' and 'T'. Its
% dimensional deficit equals 0. A one-dimensional family of segments is
% obtained by releasing 'T' and allowing it to vary during continuation.

coll_args = {@linode, @linode_DFDX, @linode_DFDP, 0, [1 0 0], 'p', 1};
prob = bvp_isol2seg(coco_prob(), '', coll_args{:}, @lin_bc);
[data uidx] = coco_get_func_data(prob, 'bvp.seg.coll', 'data', 'uidx');
prob = coco_add_pars(prob, 'duration', uidx(data.T_idx), 'T');
prob = coco_add_slot(prob, 'IP', @add_IP, [], 'bddat');
bd = coco(prob, 'run1', [], 1, 'T', [0 2]);

x0 = coco_bd_col(bd, 'X0'); % Extract column data
plot(x0(1,:),x0(2,:),'r')

coco_use_recipes_toolbox % remove the coll_v2 and bvp_v1 toolboxes from the search path
