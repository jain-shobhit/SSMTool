coco_use_recipes_toolbox varcoll_v1 bvp_v1 coll_v1 % Add bvp_v1, varcoll_v1, and coll_v1 toolboxes to search path

%% Example 10.2

% The continuation problem structure encoded below corresponds to a family
% of 151 zero functions (120 collocation conditions, 27 continuity
% conditions, and four boundary conditions) in terms of 152 continuation
% variables (150 basepoint values, one interval length, and one problem
% parameter), a single monitor function that evaluates to the problem
% parameter, and the corresponding inactive continuation parameter 'p'. Its
% dimensional deficit equals 0. A one-dimensional family of periodic orbits
% is obtained by releasing 'p' and allowing it to vary during continuation.
% The corresponding Floquet multipliers are stored with the bifurcation
% data.

% Note that the screen output during the initialization of the atlas
% algorithm differs from what is printed in Recipes for Continuation.

[t0 x0]   = ode45(@(t,x) linode(x, 1), [0 2*pi], [0; 1; 0]); % Approximate periodic orbit
coll_args = {@linode, @linode_DFDX, @linode_DFDP, t0, x0, 'p', 1};
prob = bvp_isol2seg(coco_prob(), '', coll_args{:}, @lin_bc);

prob = var_coll_add(prob, 'bvp.seg'); % Add 'varcoll' instance
prob = po_mult_add(prob, 'bvp.seg');  % Compute Floquet multipliers
bd = coco(prob, 'run', [], 1, 'p', [0.2 2]);

fm = coco_bd_col(bd, 'bvp.seg.multipliers') % Extract Floquet multipliers (including trivial at 1)

coco_use_recipes_toolbox % Remove bvp_v1, varcoll_v1, and coll_v1 toolboxes from search path
