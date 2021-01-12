coco_use_recipes_toolbox coll_v1 % add the coll_v1 toolbox to the search path

%% Section 7.3.1

% The continuation problem encoded below corresponds to a family of 98 zero
% function (80 collocation conditions, 18 continuity conditions), in terms
% of 101 continuation variables (100 basepoint values and one interval
% length), a family of four monitor functions that evaluate to the
% components of the trajectory end point at t=0, the first component of the
% trajectory end point at t=1, and the interval length, respectively, and
% four corresponding inactive continuation parameters 'y1s', 'y2s', 'y1e',
% and 'T'. Its dimensional deficit equals -1. The call to the coco
% entry-point function indicates a desired manifold dimension of 1. To this
% end, the continuation parameter 'T' and 'y1' are released and allowed to
% vary during continuation.

t0 = [0; 0.04];
x0 = [1 0; 1 0.04];

prob = coll_isol2seg(coco_prob(), '', @catenary, t0, x0, []); % Build 'coll' continuation problem
data = coco_get_func_data(prob, 'coll', 'data'); % Extract toolbox data
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.x1_idx(1); data.T_idx], ...
  {'y1s' 'y2s' 'y1e' 'T'});
coco(prob, 'coll1', [], 1, {'T' 'y1e'}, [0 1]);

% Use the following commands to graph the final discretized trajectory, see
% also figure_7_1.m.

% sol = coll_read_solution('', 'coll1', 5);
% plot(sol.t, sol.x(:,1), 'r')

% The continuation problem encoded below is identical to that above, but
% constructed from a stored solution from the previous run. The call to the
% coco entry-point function indicates a desired manifold dimension of 1. To
% this end, the continuation parameter 'y1e' and 'y2s' are released and
% allowed to vary during continuation.

prob = coll_sol2seg(coco_prob(), '', 'coll1', 5); % Reconstruct 'coll' continuation problem
data = coco_get_func_data(prob, 'coll', 'data');
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.x1_idx(1); data.T_idx], ...
  {'y1s' 'y2s' 'y1e' 'T'});
coco(prob, 'coll2', [], 1, {'y1e' 'y2s'}, [0 3]);

coco_use_recipes_toolbox % remove the coll_v1 toolbox from the search path
