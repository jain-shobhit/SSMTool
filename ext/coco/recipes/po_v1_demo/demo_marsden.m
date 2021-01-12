coco_use_recipes_toolbox coll_v1 po_v1 % add the coll_v1 and po_v1 toolboxes to the search path

%% Example 8.3

% The continuation problem encoded below corresponds to a family of 151
% zero functions (120 collocation conditions, 27 continuity conditions,
% three boundary conditions, one integral phase condition) in 153
% continuation variables (150 basepoint values, one interval length, two
% problem parameters), three monitor functions (two that evaluate to the
% problem parameters and one that evaluates to the interval length), and
% the corresponding inactive continuation parameters 'p1' and 'p2' and
% active continuation parameter 'po.period'. Its dimensional deficit equals
% 0. The call to the coco entry-point function indicates a desired manifold
% dimension of 1. To this end, the continuation parameter 'p1' is released
% and allowed to vary during continuation.

t0 = (0:2*pi/100:2*pi)';
x0 = 0.01*(cos(t0)*[1 0 -1]-sin(t0)*[0 1 0]);
p0 = [0; 6];

prob = coco_prob();
prob = po_isol2orb(prob, '', @marsden, t0, x0, {'p1' 'p2'}, p0);
coco(prob, 'po1', [], 1, {'p1' 'po.period'}, [-1 1]);

%% Example 8.4

% The continuation problem structure encoded below is identical to that
% above, but an initial solution guess is obtained by various modifications
% to a periodic orbit found in the previous run. The call to coco_xchg_pars
% constrains the interval length, while releasing the second problem
% parameter. The family of high-period periodic orbits of constant period
% approximates a family of homoclinic connections to a saddle equilibrium.

sol1 = po_read_solution('', 'po1', 6); % Periodic orbit with T=1.9553e+01
f = marsden(sol1.x', repmat(sol1.p, [1 size(sol1.x, 1)])); % Evaluate vector field at basepoints
[mn idx] = min(sqrt(sum(f.*f, 1))); % Find basepoint closest to equilibrium

scale = 25;
T  = sol1.t(end);
t0 = [sol1.t(1:idx,1) ; T*(scale-1)+sol1.t(idx+1:end,1)]; % Crank up period by factor scale
x0 = sol1.x;
p0 = sol1.p;

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', ceil(scale*10)); % Increase mesh resolution by factor scale
prob = po_isol2orb(prob, '', @marsden, t0, x0, {'p1' 'p2'}, p0);
prob = coco_xchg_pars(prob, 'p2', 'po.period'); % Exchange parameters to fix period and free p2
coco(prob, 'po2', [], 1, {'p1' 'p2' 'po.period'}, [-1 1]);

coco_use_recipes_toolbox % remove the coll_v1 and po_v1 toolboxes from the search path
